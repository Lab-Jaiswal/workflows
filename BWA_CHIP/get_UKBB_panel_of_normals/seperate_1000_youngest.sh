#!/bin/bash

#Download the 1000 youngest file list (created in get_1000_youngest_with_crams_no_predispositions.R)
dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/top_1000_without_predis_with_crams

#Get a list of all of the WES file paths in DNANexus
#Add a / before every space in the filepath
dx ls --full project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/ | xargs -I {} bash -c 'array=( {array[@]/#/{}/}"' > all_file_paths
sed -e 's/\s/\ /g' < all_file_paths > all_file_paths_fixed

#Filter out all filepaths to only list the 1000 youngest
#Add a / before every space in the filepath
#Add the projec name to every row
grep -f top_1000_without_predis_with_crams < all_file_paths_fixed > filepaths_of_1000_youngest
sed -e 's/\s/\ /g' < filepaths_of_1000_youngest > filepaths_of_1000_youngest_fixed
sed -e 's/^/project-G5B07V8JPkg740v9GjfF9PzV:/g' filepaths_of_1000_youngest_fixed > filepaths_of_1000_youngest_final

#upload the 1000 youngest filepath to dna nexus
dx upload filepaths_of_1000_youngest_final --path project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/

#dx download the 1000 youngest using the filepaths in step 5
mkdir 1000_youngest
cd 1000_youngest
cat ../filepaths_of_1000_youngest_final | xargs -I % dx download %

#upload the folder with the 1000 youngest
