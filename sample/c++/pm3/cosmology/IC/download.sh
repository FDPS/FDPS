TARGET=ic_sb128.tar
wget -O ./${TARGET} "https://v2.jmlab.jp/owncloud/index.php/s/XnzvW5XAYwfqZYQ/download?path=%2Fsb&files=${TARGET}"
tar -xf ./${TARGET}
rm ./${TARGET}

