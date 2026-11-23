param([int]$MaxSeq = 2, [int]$Wait = 3)
$root = 'D:\CAAS\Recombination\WatkinsRILsRecombinationStudy'
$fa = Join-Path $root 'Tables\Table48_candidate_proteins.faa'
$out = Join-Path $root 'Tables\Table48b_blast_results'
New-Item -ItemType Directory -Force -Path $out | Out-Null
$lines = Get-Content $fa
$seqs = @(); $id = $null; $buf = @()
foreach ($l in $lines) {
  if ($l.StartsWith('>')) {
    if ($id) { $seqs += [pscustomobject]@{ id = $id; seq = ($buf -join '') } }
    $id = $l.Substring(1); $buf = @()
  } else { $buf += $l.Trim() }
}
if ($id) { $seqs += [pscustomobject]@{ id = $id; seq = ($buf -join '') } }
"sequences: $($seqs.Count)"
$done = 0
foreach ($s in $seqs) {
  if ($done -ge $MaxSeq) { break }
  $dest = Join-Path $out ($s.id + '.txt')
  if (Test-Path $dest) { $done++; continue }
  $body = @{ CMD='Put'; PROGRAM='blastp'; DATABASE='refseq_protein';
             QUERY=$s.seq; HITLIST_SIZE='5'; EXPECT='1e-5'; FORMAT_TYPE='Text' }
  try {
    $r = Invoke-WebRequest -Uri 'https://blast.ncbi.nlm.nih.gov/Blast.cgi' -Method Post -Body $body -UseBasicParsing -TimeoutSec 120
    $rid = [regex]::Match($r.Content, 'RID = (\S+)').Groups[1].Value
    if (-not $rid) { "no RID for $($s.id)"; continue }
    $res = $null
    for ($i = 0; $i -lt 40; $i++) {
      Start-Sleep -Seconds $Wait
      $g = Invoke-WebRequest -Uri ("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$rid&FORMAT_TYPE=Text&ALIGNMENTS=5&DESCRIPTIONS=5") -UseBasicParsing -TimeoutSec 120
      if ($g.Content -match 'Sequences producing significant alignments|No hits found') { $res = $g.Content; break }
      if ($g.Content -match 'Status=WAITING') { continue }
    }
    if ($res) { Set-Content -Path $dest -Value $res -Encoding UTF8; $done++ }
  } catch { "failed $($s.id): $($_.Exception.Message)" }
}
"done: $done"