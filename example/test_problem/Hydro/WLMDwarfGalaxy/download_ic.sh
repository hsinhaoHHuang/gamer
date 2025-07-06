#!/bin/bash

LOCAL_FILENAME1="ICsForTestProbWLMDwarfGalaxy_Ng7e6"
FILE_ID1="6869527eed776d031cebfe68"

# file download
curl https://hub.yt/api/v1/item/${FILE_ID1}/download -o "${LOCAL_FILENAME1}.tar.gz"
wget https://github.com/grackle-project/grackle_data_files/raw/refs/heads/main/input/CloudyData_noUVB.h5

# file unzip
tar xzvf ${LOCAL_FILENAME1}.tar.gz
mv ${LOCAL_FILENAME1}/* ./
rmdir ${LOCAL_FILENAME1}
rm ${LOCAL_FILENAME1}.tar.gz
