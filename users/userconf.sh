#!/bin/bash

groupadd users

useradd -m -g 100 -u 1003 -s /bin/bash sangsoon
echo "sangsoon:axio" | chpasswd

useradd -m -g 100 -u 1007 -s /bin/bash yongfang
echo "yongfang:axio" | chpasswd

useradd -m -g 100 -u 1010 -s /bin/bash joseph
echo "joseph:axio" | chpasswd

useradd -m -g 100 -u 1012 -s /bin/bash beatriz
echo "beatriz:axio" | chpasswd

userdel rstudio
