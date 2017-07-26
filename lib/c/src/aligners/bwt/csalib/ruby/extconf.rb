require 'mkmf'
$objs = %w{
  ../csa.a ruby.o
}
create_makefile('csa')
