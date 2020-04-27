REM FOR %%s IN (*.nex) DO BayesTraits.V2.Dec.2014.exe %%s >> scriptout.txt
FOR %%s IN (*.nex) DO BayesTraits.V2.Dec.2014.exe %%~ns.nex %%~ns.txt <RJMCMC_commands.txt >> scriptout.txt
