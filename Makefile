DEPTH = .
SUBDIRS = misc geometry samplers restraints builders scep appscripts
TESTDIRS = catrace xcatrace protlig protrna xloop ppi ssEM

-include $(DEPTH)/makedefs

webdist:
	@echo $@
	@mkdir /tmp/webdist
	@svn export -r 103 . /tmp/webdist/rtk
	@rm -Rf /tmp/webdist/rtk/chapter
	@cp -r /tmp/webdist/rtk/html /tmp/webdist
	@cd /tmp/webdist ; tar cvzf rtk.tgz rtk
	@rm -Rf /tmp/webdist/rtk

# DO NOT DELETE
