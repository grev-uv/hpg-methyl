# This SConstruct does nothing more than load the SConscript in this dir
# The Environment() is created in the SConstruct script
# This dir can be built standalone by executing scons here, or together
# by executing scons in the parent directory
import sys

debug = int(ARGUMENTS.get('debug', '0'))

compiler = ARGUMENTS.get('compiler', 'gcc')

SConscript('SConscript', exports = ['debug', 'compiler'])
