# -*- PyInstaller input file -*-
# -*- mode: python           -*-

try:
    import pyInstaller
    pyInstaller_path = [os.path.dirname(pyinstaller.__file__)]
    print("pyInstaller_path = %r" % pyInstaller_path)
except ImportError:
    pyInstaller_path = []



import os
import sys
import shutil



python_path = os.path.dirname(sys.executable)
print('python_path ', python_path)
pkg_path = os.path.abspath(os.path.join('.', '..', 'tesis-desktop'))
print('pkg_path ', pkg_path)

pathex = pyInstaller_path + [
    python_path,
    os.path.join(python_path, 'Lib'),
    os.path.join(python_path, 'Lib', 'site-packages'),
]

mkl_dlls = [
    os.path.join(python_path, 'evns', 'py27', 'Library', 'bin', 'mkl_def3.dll')
]

has_mkl_dlls = False
if mkl_dlls:
    mkl_dll_base = os.path.basename(mkl_dlls[0])
    #assert os.path.exists(mkl_dll), '%s doesnt exist' % mkl_dll
    has_mkl_dlls = os.path.exists(mkl_dll_base)

binaries = []
if sys.platform == 'win32':
    binaries = [
        ('msvcp100.dll', 'C:\\Windows\\System32\\msvcp100.dll', 'BINARY'),
        ('msvcr100.dll', 'C:\\Windows\\System32\\msvcr100.dll', 'BINARY'),
    ]
    if has_mkl_dlls:
        for mkl_dll in mkl_dlls:
            binaries.append(
                (mkl_dll_base, mkl_dll, 'BINARY')
            )


hiddenimports = ['Bio.SearchIO.BlastIO']

excludes=['FixTk', 'tcl', 'tk', '_tkinter', 'tkinter', 'Tkinter']


block_cipher = None
added_files=[('./project/Blastdb','Blastdb'),
('./project/BlastResult','Categories'),
('./project/Categories','BlastResult'),
('./project/Databases','Databases'),
('./project/DbAmbigua','DbAmbigua'),
('./project/FinalResult','FinalResult'),
('./project/Test','Test'),
('./project/tmp.fa','.'),
('./project/View', 'View')
  ]


a = Analysis(['tesis.py'],
             datas=added_files,
             hookspath=[],
             runtime_hooks=[],
             pathex=pathex,
             excludes=excludes,
             hiddenimports=hiddenimports,
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries + binaries,
          a.zipfiles,
          a.datas,
          name='tesis',
          debug=True,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='tesis')
