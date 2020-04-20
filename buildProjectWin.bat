rmdir /s /q build 
rmdir /s /q bin 

mkdir build
mkdir build

cd build

cmake -G "MinGW Makefiles" ..

REM make all
mingw32-make -j4
cd ..
