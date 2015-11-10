#!/bin/sh
# path to comparison program
MSSIM="/home/student/tabkhi/project/bge/ms_ssim/ms_ssim_bmp"
# path to reference output
GOLD="../sequ-result"

   fgtotal=0.0
   bgtotal=0.0
   for i in {150..197} 
     do 
       FG_SSIM=`$MSSIM fg0$i.bmp $GOLD/fg0$i.bmp`
       fgtemp=$(echo "scale=2; ${FG_SSIM} + ${fgtotal}" | bc)
       fgtotal=$fgtemp
	echo $FG_SSIM

       BG_SSIM=`$MSSIM bg0$i.bmp $GOLD/bg0$i.bmp`
       bgtemp=$(echo "scale=2; ${BG_SSIM} + ${bgtotal}" | bc)
       bgtotal=$bgtemp
     done

   fgfinal=$(echo "scale=2; ${fgtotal} / 48.00" | bc)
   bgfinal=$(echo "scale=2; ${bgtotal} / 48.00" | bc)

   echo "fg: $fgfinal    bg: $bgfinal" >> ./res.txt
