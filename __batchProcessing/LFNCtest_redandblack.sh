#!/bin/bash

MAINDIR=$( dirname $( cd "$( dirname $0 )" && pwd ) );
TMCDIR=${MAINDIR}/mpeg-pcc-tmc2-netTest

## Input parameters
SRCDIR=${MAINDIR}/pointCloudCompress_testSequence
EXTERNAL=${MAINDIR}/external;
EXEDIR=${TMCDIR}/bin;
CFGDIR=${TMCDIR}/cfg

SEQ=24;       # in [22;26]
FRAMECOUNT=32;
MANNER="MLPnoCBFOM_03_06";
COND="C2AI";       # in [C2AI, C2RA, CWAI, CWRA]
THREAD=1;
METRICRESOLUTION=1023;
if [ $SEQ == 27 -o $SEQ == 28 ]; then METRICRESOLUTION=2047; fi
if [ $SEQ == 29 ]; then METRICRESOLUTION=4095; fi

## Set external tool paths
NDATAPATH="set_value";
SEQNAME="set_value";
OUTPUT="set_value";
ENCODER=${EXEDIR}/${MANNER}/PccAppEncoder.exe;
HDRCONVERT=${EXTERNAL}/HDRTools-v0.18/bin/HDRConvert.exe;

## Parameters and pathes check
if [ ! -f $ENCODER    ] ; then echo "Can't find PccAppEncoder, please set.     ($ENCODER )";   exit -1; fi
if [ ! -f $HDRCONVERT ] ; then echo "Can't find HdrConvert, please set.        ($HDRCONVERT)"; exit -1; fi

## Set Configuration based on sequence, condition and rate
if [ $COND == "C2AI" -o $COND == "C2RA" ] 
then
  case $SEQ in
      22) CFGSEQUENCE="sequence/queen.cfg";;
      23) CFGSEQUENCE="sequence/loot_vox10.cfg";;
      24) CFGSEQUENCE="sequence/redandblack_vox10.cfg";;
      25) CFGSEQUENCE="sequence/soldier_vox10.cfg";;
      26) CFGSEQUENCE="sequence/longdress_vox10.cfg";;
      27) CFGSEQUENCE="sequence/basketball_player_vox11.cfg";;
      28) CFGSEQUENCE="sequence/dancer_vox11.cfg";;
      29) CFGSEQUENCE="sequence/Thaidancer_vox12.cfg";;
      *) echo "sequence not correct ($SEQ)";   exit -1;;
  esac
else
   case $SEQ in
      22) CFGSEQUENCE="sequence/queen-lossless.cfg";;
      23) CFGSEQUENCE="sequence/loot_vox10-lossless.cfg";;
      24) CFGSEQUENCE="sequence/redandblack_vox10-lossless.cfg";;
      25) CFGSEQUENCE="sequence/soldier_vox10-lossless.cfg";;
      26) CFGSEQUENCE="sequence/longdress_vox10-lossless.cfg";;
      27) CFGSEQUENCE="sequence/basketball_player_vox11-lossless.cfg";;
      28) CFGSEQUENCE="sequence/dancer_vox11-lossless.cfg";;
      29) CFGSEQUENCE="sequence/Thaidancer_vox12-lossless.cfg";;
      *) echo "sequence not correct ($SEQ)";   exit -1;;
  esac 
  CFGRATE="rate/ctc-r5.cfg"
  BIN=S${SEQ}${COND}_F${FRAMECOUNT}.bin
fi

case $COND in
  CWAI) CFGCOMMON="common/ctc-common-lossless-geometry-texture.cfg";;
  CWLD) CFGCOMMON="common/ctc-common-lossless-geometry-texture.cfg";; 
  C2AI) CFGCOMMON="common/ctc-common.cfg";;                           
  C2RA) CFGCOMMON="common/ctc-common.cfg";;           
  *) echo "Condition not correct ($COND)";   exit -1;;
esac

case $COND in
  CWAI) CFGCONDITION="condition/ctc-all-intra-lossless-geometry-texture.cfg";;
  CWLD) CFGCONDITION="condition/ctc-low-delay-lossless-geometry-texture.cfg";;
  C2AI) CFGCONDITION="condition/ctc-all-intra.cfg";;
  C2RA) CFGCONDITION="condition/ctc-random-access.cfg";;  
  *) echo "Condition not correct ($COND)";   exit -1;;
esac

case $SEQ in
  22) NDATAPATH="queen/queen_n/frame_%04d_n.ply";SEQNAME="queen";;
  23) NDATAPATH="loot/loot_n/loot_vox10_%04d_n.ply";SEQNAME="loot";;
  24) NDATAPATH="redandblack/redandblack_n/redandblack_vox10_%04d_n.ply";SEQNAME="redandblack";;
  25) NDATAPATH="soldier/soldier_n/soldier_vox10_%04d_n.ply";SEQNAME="soldier";;
  26) NDATAPATH="longdress/longdress_n/longdress_vox10_%04d_n.ply";SEQNAME="longdress";;
  27) NDATAPATH="basketball_player/basketball_player_vox11_%08d.ply";SEQNAME="basketball_player";;
  28) NDATAPATH="dancer/dancer_vox11_%08d.ply";SEQNAME="dancer";;
  29) NDATAPATH="Thaidancer/Thaidancer_n/Thaidancer_viewdep_vox12_%06d_n.ply";SEQNAME="Thaidancer";;
  *) echo "sequence not correct ($SEQ)";   exit -1;;
esac 

## Encoder 
for RATE in {1..5}
do
  case $RATE in
      5) CFGRATE="rate/ctc-r5.cfg";OUTPUT=${MAINDIR}"/__output/"${MANNER}"/"${SEQNAME}${FRAMECOUNT}_${COND}"/r5/";;
      4) CFGRATE="rate/ctc-r4.cfg";OUTPUT=${MAINDIR}"/__output/"${MANNER}"/"${SEQNAME}${FRAMECOUNT}_${COND}"/r4/";;
      3) CFGRATE="rate/ctc-r3.cfg";OUTPUT=${MAINDIR}"/__output/"${MANNER}"/"${SEQNAME}${FRAMECOUNT}_${COND}"/r3/";;
      2) CFGRATE="rate/ctc-r2.cfg";OUTPUT=${MAINDIR}"/__output/"${MANNER}"/"${SEQNAME}${FRAMECOUNT}_${COND}"/r2/";;
      1) CFGRATE="rate/ctc-r1.cfg";OUTPUT=${MAINDIR}"/__output/"${MANNER}"/"${SEQNAME}${FRAMECOUNT}_${COND}"/r1/";;
      *) echo "rate not correct ($RATE)";   exit -1;;
  esac
  BIN=${OUTPUT}S${SEQ}${COND}R0${RATE}_F${FRAMECOUNT}.bin
  
  $ENCODER \
    --config=${CFGDIR}/${CFGCOMMON} \
    --config=${CFGDIR}/${CFGSEQUENCE} \
    --config=${CFGDIR}/${CFGCONDITION} \
    --config=${CFGDIR}/${CFGRATE} \
    --configurationFolder=${CFGDIR}/ \
    --absoluteD1=1 \
    --absoluteT1=1 \
    --keepIntermediateFiles=1 \
    --uncompressedDataFolder=${SRCDIR}/ \
    --frameCount=$FRAMECOUNT \
    --colorSpaceConversionPath=$HDRCONVERT \
    --reconstructedDataPath=${BIN%.???}_rec_%04d.ply \
    --resolution=$METRICRESOLUTION \
    --compressedStreamPath=$BIN \
    --normalDataPath=${SRCDIR}/$NDATAPATH | tee -a ${MAINDIR}/__statisticData/${MANNER}${FRAMECOUNT}_${SEQNAME}_${COND}.txt

done