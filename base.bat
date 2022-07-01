@echo off

echo case%1-------------------------------------------

copy Data\Unst_case%1.dat         Data\Unst.dat
copy Data\Vert_Bottom_case%1.dat  Data\Vert_Bottom.dat

copy Data\Wind_%2.dat             Data\Wind.dat
copy Data\WL_%2.dat               Data\WL.dat

copy 3D_Model_Simple_%3.exe       3D_Model_Simple.exe

call 3D_Model_Simple.exe

del Data\Unst.dat
del Data\Vert_Bottom.dat

del Data\Wind.dat

del 3D_Model_Simple.exe

