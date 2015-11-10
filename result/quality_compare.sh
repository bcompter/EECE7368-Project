#!/bin/sh
# path to comparison program
MSSIM="/home/student/tabkhi/project/bge/ms_ssim/ms_ssim_bmp"
# path to reference output
GOLD="./"

   fgtotal=0.0
   bgtotal=0.0
   for i in {1..1} 
     do 
#FG_SSIM=`$MSSIM 9-window-$i.bmp $GOLD/11-window-$i.bmp`
       FG_SSIM=`$MSSIM 3_window.bmp $GOLD/11_window.bmp`
       echo $FG_SSIM
       fgtemp=$(echo "scale=2; ${FG_SSIM} + ${fgtotal}" | bc)
       fgtotal=$fgtemp
	echo $FG_SSIM

#       BG_SSIM=`$MSSIM bg0$i.bmp $GOLD/bg0$i.bmp`
#       bgtemp=$(echo "scale=2; ${BG_SSIM} + ${bgtotal}" | bc)
#       bgtotal=$bgtemp
     done

   fgfinal=$(echo "scale=2; ${fgtotal} / 1.00" | bc)
#bgfinal=$(echo "scale=2; ${bgtotal} / 48.00" | bc)

#   echo "fg: $fgfinal    bg: $bgfinal" >> ./res.txt
   echo "fg: $fgfinal " 
