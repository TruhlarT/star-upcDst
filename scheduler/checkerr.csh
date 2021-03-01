#!/bin/csh


ls -l /star/data01/pwg_tasks/upc02/Parta/logs/*.err | awk '{print $5, $9}' | grep -v "213 /s"