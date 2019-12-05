rmdir /s /q build 
rmdir /s /q bin 

mkdir build
mkdir build

cd build

cmake -G "MinGW Makefiles" ..

make all

cd ..