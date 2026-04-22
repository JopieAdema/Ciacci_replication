@echo off
rem Run both replication scripts in order.
rem Logs land in C:/Users/c4041171/replication.log and
rem              C:/Users/c4041171/Dropbox/.../output/replication_fuzzydid.log.

set RSCRIPT="C:/Users/c4041171/AppData/Local/Programs/R/R-4.4.3/bin/Rscript.exe"
set HELPER="C:/Users/c4041171/run_stata.R"
set ROOT=C:/Users/c4041171/Dropbox/Ciacci replication/replication_220426

echo === Running replication.do ===
%RSCRIPT% %HELPER% "%ROOT%/code/replication.do"

echo === Running replication_fuzzydid.do ===
%RSCRIPT% %HELPER% "%ROOT%/code/replication_fuzzydid.do"

echo === All done. Check %ROOT%/output/ ===
