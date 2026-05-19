#!/usr/bin/env pwsh
# ! --- response_spectrum_mystran_add begin --- !
param(
  [Parameter(Mandatory = $true)]
  [string]$Bdf,

  [ValidateSet("FAST", "CHECK", "FEMAP")]
  [string]$Mode = "FEMAP",

  [string]$MystranExe = "E:\mystran17\mystran\Binaries\mystran.exe",
  [string]$PythonExe = "C:\Users\arypr\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe",

  [switch]$WriteNeu,
  [string]$GeometryNeu,
  [string]$ModalNeu,
  [string]$ModalF06,
  [string]$RsxNeu,
  [string]$RsyNeu,
  [string]$OutNeu
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if (!(Test-Path $MystranExe)) { throw "MYSTRAN executable not found: $MystranExe" }
if (!(Test-Path $Bdf)) { throw "BDF not found: $Bdf" }

Write-Host "Mode      : $Mode"
Write-Host "Solver    : $MystranExe"
Write-Host "Input BDF : $Bdf"

switch ($Mode) {
  "FAST"  { Write-Host "Profile   : PLOT-heavy (minimal F06 text)"; break }
  "CHECK" { Write-Host "Profile   : PRINT+PLOT (verbose validation)"; break }
  "FEMAP" { Write-Host "Profile   : PLOT + optional NEU post-step"; break }
}

$sw = [System.Diagnostics.Stopwatch]::StartNew()
& $MystranExe $Bdf
$solverExit = $LASTEXITCODE
$sw.Stop()
Write-Host ("Solve elapsed (s): {0:N3}" -f $sw.Elapsed.TotalSeconds)

if ($solverExit -ne 0) {
  throw "MYSTRAN returned non-zero exit code: $solverExit"
}

$baseNoExt = [System.IO.Path]::ChangeExtension($Bdf, $null)
$f06Path = "$baseNoExt.F06"
$errPath = "$baseNoExt.ERR"

if (Test-Path $f06Path) {
  $f06Tail = (Get-Content $f06Path -Tail 40) -join "`n"
  if ($f06Tail -notmatch 'terminated normally') {
    throw "MYSTRAN did not report normal termination in F06 tail: $f06Path"
  }
}

if (Test-Path $errPath) {
  $errTail = (Get-Content $errPath -Tail 120) -join "`n"
  if (($errTail -match '\*ERROR') -or ($errTail -match 'FATAL')) {
    throw "ERR tail indicates error/fatal. Check: $errPath"
  }
}

if (-not $WriteNeu) {
  Write-Host "NEU step  : skipped (use -WriteNeu to generate FEMAP NEU)"
  exit 0
}

if (!(Test-Path $PythonExe)) { throw "Python runtime not found: $PythonExe" }

$missing = @()
foreach ($p in @($GeometryNeu, $ModalNeu, $ModalF06, $RsxNeu, $RsyNeu)) {
  if ([string]::IsNullOrWhiteSpace($p) -or !(Test-Path $p)) { $missing += $p }
}
if ($missing.Count -gt 0) {
  throw "Missing NEU writer inputs. Provide valid: -GeometryNeu -ModalNeu -ModalF06 -RsxNeu -RsyNeu"
}
if ([string]::IsNullOrWhiteSpace($OutNeu)) {
  throw "Missing -OutNeu output path."
}

$writer = "E:\mystran17\codex_mod\response_spectrum_mystran_add\scripts\femap_rs_neu_writer.py"
if (!(Test-Path $writer)) { throw "NEU writer not found: $writer" }

Write-Host "NEU step  : running femap_rs_neu_writer.py"
& $PythonExe $writer `
  --geometry $GeometryNeu `
  --modal-neu $ModalNeu `
  --modal-f06 $ModalF06 `
  --rsx-neu $RsxNeu `
  --rsy-neu $RsyNeu `
  --out $OutNeu

Write-Host "NEU output: $OutNeu"
# ! --- response_spectrum_mystran_add end --- !
