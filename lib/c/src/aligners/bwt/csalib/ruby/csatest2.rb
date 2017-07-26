require 'csa'
require 'benchmark'

def traverse(csa, range, th, key)
#  print("traverse [", range[0], ",", range[1], "] key = ", key, "\n")
  f = 0
  csa.child_l(range).each {|c, range2|
    num = range2[1] - range2[0] + 1
    if num >= th then
      newkey = "%c" % c;
      newkey += key
      traverse(csa, range2, th, newkey)
      f = 1
    end
  }
  if f == 0 then
    t = csa.child_r(key.length, range)
    if t.size != 1 then
      print("[", range[0], ",", range[1], "] key = ", key, "\n")
    end
  end
end

def traverse2(csa, range, th, key, depth)
#  print("traverse2 [", range[0], ",", range[1], "] key = ", key, "\n")
  f = 0
  csa.child_r(depth,range).each {|c, range2|
    num = range2[1] - range2[0] + 1
    if num >= th then
      newkey = "%c" % c;
      newkey = key + newkey;
      traverse2(csa, range2, th, newkey, depth+1)
      f = 1
    end
  }
  if f == 0 then
    print("[", range[0], ",", range[1], "] key = ", key, "\n")
  end
end


csa = CSA.new(ARGV[0],ARGV[1])

n = csa.getn()
th = n / 10000

puts Benchmark.measure{
traverse(csa, [0, n], th, "")
#traverse2(csa, [0, n], th, "", 0)
}
