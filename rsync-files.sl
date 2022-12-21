#!/bin/bash
rsync -v --exclude-from=rsync-exclude-list.txt -r --delete . mdm95@mahuika:/nesi/project/uc03144/packages/o3seespy
# rsync -v --exclude-from=rsync-exclude-list.txt -r --delete . mdm95@mahuika:packages/o3seespy