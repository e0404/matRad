#!/bin/bash
$1/topas_extensions
wget -O $1/topas_extensions/extensions.zip https://github.com/topasmc/extensions/archive/master.zip
unzip $1/topas_extensions/extensions.zip -d $1/topas_extensions
mv $1/topas_extensions/extensions-master/* $1/topas_extensions/
rm -r $1/topas_extensions/extensions-master
rm $1/topas_extensions/extensions.zip
