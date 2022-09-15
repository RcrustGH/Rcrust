# Makes backup of old .RData file
if [ -f .RData ]; then
	mv .RData .RData.backup
fi

# Copies the Rcrust RData file to the default to be loaded, and then
# executes the Rcrust program GUI
cp Rcrust.RData .RData
# Creates temporary file for quit command to close interactive terminal at end
echo "q()" > .tmp_cmd
R --restore --no-restore-history --q --interactive --no-save < .tmp_cmd
rm .tmp_cmd

# Restores old .RData file (if it existed) and cleans up
rm .RData
if [ -f .RData.backup ]; then
	mv .RData.backup .RData 
fi
