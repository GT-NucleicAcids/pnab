setlocal EnableDelayedExpansion

CALL "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" x86_amd64
CALL set PATH=%PATH%;C:\Users\alenaizan\Miniconda3\condabin\
CALL conda activate pnab
IF NOT EXIST build\ CALL cmake -G "NMake Makefiles" -Bbuild -DCMAKE_BUILD_TYPE:STRING=Release -DOPENBABEL_DIR="%CONDA_PREFIX%"
CALL cd build
CALL nmake
CALL XCOPY /Y bind* ..\pnab
CALL cd ..
CALL set BABEL_DATADIR=C:\Users\alenaizan\Miniconda3\envs\pnab\share\openbabel
CALL set SP_DIR="%CONDA_PREFIX%\Lib\site-packages"
CALL XCOPY /E /I /Y pnab "%SP_DIR%\pnab"
CALL XCOPY /E /I /Y tests "%SP_DIR%\pnab\tests"
CALL dir "%SP_DIR%\pnab\"
CALL cd tests
CALL pytest -s