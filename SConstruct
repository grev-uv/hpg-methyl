
import os

# Initialize the environment with path variables, CFLAGS, and so on
hpglib_path = '#lib/c/'
third_party_path = '#/lib/third_party'
third_party_samtools_path = '#/lib/third_party/samtools'
third_party_hts_path = '#/lib/third_party/htslib'

vars = Variables('buildvars.py')

compiler = ARGUMENTS.get('compiler', 'gcc')
build_tools = [ 'default', 'packaging' ]

env = Environment(tools = build_tools,
  CC = compiler,
  CFLAGS = '-Wall -std=c99 -D_GNU_SOURCE -fopenmp -D_REENTRANT -msse4.2',
  CPPPATH = ['#', '#src', hpglib_path + 'src', third_party_path, third_party_samtools_path, third_party_hts_path, '/usr/include', '/usr/local/include', '/usr/include/libxml2'],
  LIBPATH = [hpglib_path + 'build', third_party_hts_path, third_party_samtools_path, '/usr/lib', '/usr/local/lib'],
  LIBS = ['curl', 'dl', 'gsl', 'gslcblas', 'm', 'xml2', 'z'],
  LINKFLAGS = ['-fopenmp'])


if int(ARGUMENTS.get('debug', '0')) == 1:
  env['CFLAGS'] += ' -O0 -g'
else:
  env['CFLAGS'] += ' -O3 -g'

env['objects'] = []


# Targets
SConscript(['lib/SConstruct'])

env.Program('#bin/hpg-methyl',
  source = [
    Glob('src/*.c'),
  	Glob('src/build-index/*.c'),
  	Glob('src/bs/*.c'),
    "%s/build/libhpg.a" % hpglib_path,
    "%s/libbam.a" % third_party_samtools_path,
    "%s/libhts.a" % third_party_hts_path
  ]
)
