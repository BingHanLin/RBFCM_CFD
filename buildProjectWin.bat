rmdir /s /q build 
rmdir /s /q bin 

mkdir build
mkdir build

cd build

cmake -DCATCH2TEST=ON -G "MinGW Makefiles" ..

REM make all
mingw32-make -j4

REM run unit test
cd ..
cd bin
MYAPP_Test.exe
cd ..
