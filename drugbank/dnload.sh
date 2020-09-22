#!/bin/bash

# download the latest drugbank data
# need to specify ${user} for username and ${pswd} for password

# full database (xml file)
curl -Lfv -o fulldb.zip -u ${user}:${pswd} https://www.drugbank.ca/releases/latest/downloads/all-full-database
# xsd file: schema definition for the xml file
curl -L -o drugbank.xsd https://www.drugbank.ca/docs/drugbank.xsd
# targets, enzymes, carries and transporters for all drugs
curl -Lfv -o targets_all.zip -u ${user}:${pswd} https://www.drugbank.ca/releases/latest/downloads/target-all-polypeptide-ids
curl -Lfv -o enzymes_all.zip -u ${user}:${pswd} https://www.drugbank.ca/releases/latest/downloads/enzyme-all-polypeptide-ids
curl -Lfv -o carriers_all.zip -u ${user}:${pswd} https://www.drugbank.ca/releases/latest/downloads/carrier-all-polypeptide-ids
curl -Lfv -o transporters_all.zip -u ${user}:${pswd} https://www.drugbank.ca/releases/latest/downloads/transporter-all-polypeptide-ids

# unzip and rename
unzip fulldb.zip
mv "full database.xml" drugbank.xml
unzip targets_all.zip -d targets
unzip enzymes_all.zip -d enzymes
unzip carriers_all.zip -d carriers
unzip transporters_all.zip -d transporters

# note: for the 5.0.4 release version of the full database, the beginning of the root node is missing (but the closing, </drugbank> is there at the end of the file), so I manually add it:
#sed -i '1i <drugbank>\n' drugbank.xml
## initially I added the following, based on the data file of a previous release. But strangely with it the xml2 package in R doesn't work...
#sed -i '1i <?xml version="1.0" encoding="UTF-8"?>\n<drugbank xmlns="http://www.drugbank.ca" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.drugbank.ca http://www.drugbank.ca/docs/drugbank.xsd" version="5.0" exported-on="2017-01-01">\n' drugbank.xml
