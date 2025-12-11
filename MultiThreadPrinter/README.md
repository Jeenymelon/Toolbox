```shell
git clone https://gitee.com/anold/log4cplus.git -b 2.0.x
cd ./log4cplus
git submodule update --init --recursive
./configure --prefix=/usr/local --enable-release-version=yes --enable-thread=yes
make -j6 && sudo make install
```


