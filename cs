# Script to run ChemSol
# run command: cs molecule_name >& molecule_name.log &
  setenv SOLVINP $1.cs
  cs21.exe 
  echo ------ accounting info --------
  date
  time
  exit
