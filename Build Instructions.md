# Build Instructions
1. Clone the github repository **recursively** to get a copy of the submodules with the source: `git clone --recursive -j8 git@github.com:KingfisherA/SCESS.git`
2. Make both `Scripts\Build.sh` and `Libraries\Build_Libraries.sh` executable i.e. `chmod +x .\Scripts\Build.sh` and `chmod +x .\Libraries\Build_Libraries.sh`
3. Run `.\Scripts\Build.sh`
4. Make a cup of tea while everything compiles (it actually takes FOREVER; this is LibXC and not my terrible code I promise)
5. Be happy when you have all your binaries (yay! well done :) ) in the `.\bin` folder
6. You should be able to run `.\bin\linux\Release-x64\ViewerSCESS` straight away and the assets directory should be automatically found 
