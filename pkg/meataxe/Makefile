
# GAP_PATH=..
# PKG_PATH=${GAP_PATH}/pkg
# SRCDISTFILE=pkg

GAP_PATH=../..
PKG_PATH=..
SRCDISTFILE=meataxe

# for GAP distribution
src_dist:
	@(cmp ${PKG_PATH}/meataxe/doc/mtx.tex ${GAP_PATH}/doc/mtx.tex \
         || echo \
         "*** WARNING: current 'mtx.tex' and 'doc/mtx.tex' differ ***")
	@zoo ah ${SRCDISTFILE}.zoo \
                ${PKG_PATH}/meataxe/Makefile \
                ${PKG_PATH}/meataxe/README \
                ${PKG_PATH}/meataxe/COPYING \
                ${PKG_PATH}/meataxe/bin/Makefile \
                ${PKG_PATH}/meataxe/bin/Makefile.bc \
                `find ${PKG_PATH}/meataxe/bin/proggy -print` \
                `find ${PKG_PATH}/meataxe/doc -name "*.tex" -print` \
                `find ${PKG_PATH}/meataxe/gap -name "*.g" -print` \
                ${PKG_PATH}/meataxe/init.g \
                `find ${PKG_PATH}/meataxe/lib -print` \
                `find ${PKG_PATH}/meataxe/src -name "*.h" -print` \
                `find ${PKG_PATH}/meataxe/src -name "*.c" -print` \
                `find ${PKG_PATH}/meataxe/tests -print` 
	@zoo PE ${SRCDISTFILE}.zoo

