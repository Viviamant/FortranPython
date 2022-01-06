PYTHON=python3.8
FORTRAN=gfortran

if [ -d "build" ]; then
	rm -r build
fi

echo 'Entering "fort_src" directory'
cd fort_src

echo 'Compiling "mylib.f" with gfrotran'
$FORTRAN -fPIC -shared mylib.f -o mylib.so

echo 'Copying shared object in "/usr/local/lib"'
cp mylib.so /usr/local/lib/

cd ..

$PYTHON setup.py build
$PYTHON setup.py install
