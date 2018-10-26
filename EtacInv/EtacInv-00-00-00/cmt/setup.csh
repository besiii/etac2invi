# echo "Setting TestRelease TestRelease-00-00-80 in /afs/.ihep.ac.cn/bes3/offline/Boss/6.6.4.p03"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/contrib/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=TestRelease -version=TestRelease-00-00-80 -path=/afs/.ihep.ac.cn/bes3/offline/Boss/6.6.4.p03  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

