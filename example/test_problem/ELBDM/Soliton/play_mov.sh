for ((n=0;n<=$1;n=n+$2))
do
   echo $n
   if [ $n -lt 10 ];then
      display 'Data_00000'$n'_Projection_z_density_density.png'
      display 'Data_00000'$n'_Slice_z_density.png'
   else
      display 'Data_0000'$n'_Projection_z_density_density.png'
      display 'Data_0000'$n'_Slice_z_density.png'
   fi
done
animate -delay 100 -resize 60% *.png
