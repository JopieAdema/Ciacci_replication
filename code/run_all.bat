@echo off
rem Run both replication scripts in order, via the run_stata.R helper.
rem
rem >>> EDIT THE THREE PATHS BELOW for your machine. <<<
rem   RSCRIPT : full path to Rscript.exe
rem   HELPER  : full path to the run_stata.R wrapper
rem   ROOT    : absolute path to the replication package root
rem
rem Logs land in %ROOT%/replication.log and %ROOT%/replication_fuzzydid.log.

set RSCRIPT="C:/Users/c4041171/AppData/Local/Programs/R/R-4.4.3/bin/Rscript.exe"
set HELPER="C:/Users/c4041171/run_stata.R"
set ROOT=C:/Users/c4041171/Dropbox/Ciacci replication/replication_220426

echo === Running replication.do ===
%RSCRIPT% %HELPER% "%ROOT%/code/replication.do"

echo === Running replication_fuzzydid.do ===
%RSCRIPT% %HELPER% "%ROOT%/code/replication_fuzzydid.do"

echo === All done. Check %ROOT%/output/ ===
