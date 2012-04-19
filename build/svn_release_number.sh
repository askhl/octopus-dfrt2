#!/bin/sh
# 
# $Id: svn_release_number.sh 7117 2010-12-01 06:18:49Z dstrubbe $

cd `dirname "$0"`/..
if [ -x "$(which svn)" ] && svn info > /dev/null 2>&1 ; then
	svn info | grep Revision | awk -F: '{print $2}' | tr -d [:space:]
else
	find . -type f ! -name ChangeLog ! -name \*.svn\* \
	    ! -name \*.o ! -name \*.a ! -name \*.so \
            -exec grep '$Id:' \{\} \; \
            | grep '$Id: ' \
            | tr -d \#\!\* | awk '{print $3,"["$2,$4"]"}' \
            | grep -v 'qw(' | sort -n | tail -1

fi
