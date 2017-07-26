require 'csa'
require 'benchmark'

def dp(p1,p2, prn)
#  print("dp p1=#{p1} p2=#{p2}\n")
  n1 = p1.length
  n2 = p2.length
  d = Array.new(n2+1) {
    Array.new(n1+1,0)
  }
  d[0][0] = 0
  for i in 1..n2 do d[i][0] = i end
  for i in 1..n1 do d[0][i] = i end

  for j in 1..n2 do
    for i in 1..n1 do
      d[j][i] = [d[j-1][i-1] + ((p2[j-1] != p1[i-1])?1:0),
                 d[j-1][i]+1,
                 d[j][i-1]+1].min
    end
  end
  e = d[n2][n1]
#  p d
#  print("error = #{e}\n")
  e1 = []
  e2 = []
  i = n1
  j = n2
  while i>0 || j>0 do
    if d[j][i] == (d[j-1][i-1] + ((p2[j-1] != p1[i-1])?1:0)) then
      if p2[j-1] == p1[i-1] then
        e1 << 0
        e2 << 0
      else
        e1 << 2
        e2 << 2
      end
      i -= 1
      j -= 1
    elsif d[j][i] == d[j-1][i]+1 then
      e1 << 0
      e2 << 1
      j -= 1
    else
      e1 << 1
      e2 << 0
      i -= 1
    end
  end
#  p e1
#  p e2
if prn == 1 then
  e1.reverse!
  e2.reverse!
  i = 0
  e2.each {|x|
    if x == 2 then
      print("[","%c" % p1[i],"]")
      i += 1
    elsif x == 1 then
      print("-")
    else
      print("%c" % p1[i])
      i += 1
    end
  }
  puts
  j = 0
  e1.each {|x|
    if x == 2 then
      print("[","%c" % p2[j],"]")
      j += 1
    elsif x == 1 then
      print("-")
    else
      print("%c" % p2[j])
      j += 1
    end
  }
  puts
end
  return e
end


def approx_sub2(csa, query, error, d, range, depth, dofs, lr)
  depth += 1
  q = query.length
  d2 = Array.new(q+1)
  d2[0] = d[0]+1

  ans = []

  csa.child(dofs+depth-1,range,lr) do |c, range2|
    $nnodes += 1
    i = 1
    while i <= q do
      if lr == 0 then
        cc = query[q-i]
      else
        cc = query[i-1]
      end
      d2[i] = [d[i-1] + ((c != cc)?1:0),
               d[i]+1,
               d2[i-1]+1].min
      i += 1
    end
    if d2[q] <= error then
      ans << [d2[q],dofs+depth,range2]
    end
    if d2.min <= error then
      ans2 = approx_sub2(csa, query, error, d2, range2, depth, dofs, lr)
      if ans2 != [] then
        ans.concat(ans2)
      end
    end
  end
  return ans
end

def merge(csa, k, list)
  list.sort!{|a,b| [a[2],a[0]] <=> [b[2],b[0]]}
  list2 = []
  for i in 0..(list.length-1) do
#    print("[#{list[i][0]},#{list[i][1]},[#{list[i][2][0]},#{list[i][2][1]}]]\n")
    if i==0 || list[i][2] != list[i-1][2] then
      list2 << list[i]
    end
  end
  h = Hash.new
  list2.each{|e, d, range|
    h[range[0]] = [e, d, range]
  }

  list3 = []
  list2.each{|e, d, range|
#    print("e = #{e} range = [#{range[0]},#{range[1]}]\n")
    i = range[0]
    m = [e, d, range]
    for j in 0..k do
#      print("h[#{i}] = ")
#      p h[i]
#      print("\n")
      if h[i] == nil then
        break
      end
      if h[i][0] < e then
        m = h[i]
      end
      i = csa.lf(i)
    end
    i = range[0]
    for j in 0..k do
      if h[i] == nil then
        break
      end
      if h[i][0] < e then
        m = h[i]
      end
      i = csa.psi(i)
    end
    if m[0] == e then
      list3 << [e, d, range]
    end
  }
  list3.sort!{|a,b| [a[0],a[2]] <=> [b[0],b[2]]}
  return list3
end


def approx(csa, query, error)
  d = Array.new(query.length+1){|i| i}
  range = [0,csa.length()]
  ans = approx_sub2(csa, query, error, d, range, 0, 0, 0)

#  ans.sort!{|a,b| [a[0],a[2],a[1]] <=> [b[0],b[2],b[1]]}
#  ans.uniq!
  ans = merge(csa, error, ans)

  return ans
end

def approx2(csa, query, error)

  len = query.length

  l1 = len/2
  e1 = error/2
  q1 = query[0..l1-1]

  l2 = len - l1
  e2 = error - e1
  q2 = query[l1..len-1]

  print("l1 = #{l1} e1 = #{e1} l2 = #{l2} e2 = #{e2}\n")

  list = []

  d = Array.new(l1+1){|i| i}
  d2 = Array.new(l2+1)
  range = [0,csa.length()]
  list1 = approx_sub2(csa, q1, e1, d, range, 0, 0, 0)
  list1.each {|e, dofs, range2|
    for i in 0..l2 do d2[i] = e+i end
    tmp = approx_sub2(csa, q2, error, d2, range2, 0, dofs, 1)
    if tmp != [] then list.concat(tmp) end
  }

  d = Array.new(l2+1){|i| i}
  d1 = Array.new(l1+1)
  range = [0,csa.length()]
  list2 = approx_sub2(csa, q2, e2, d, range, 0, 0, 0)
  list2.each {|e, dofs, range2|
    for i in 0..l1 do d1[i] = e+i end
    tmp = approx_sub2(csa, q1, error, d1, range2, 0, dofs, 0)
    if tmp != [] then list.concat(tmp) end
  }

  ans = []
  list.each {|e, len, range2|
    ee = dp(query,csa.substring(range2[0], len),0)
    ans << [ee, len, range2]
  }
#  ans.sort!{|a,b| [a[0],a[2],a[1]] <=> [b[0],b[2],b[1]]}
#  ans.uniq!
  ans = merge(csa, error, ans)

  return ans
end

def approx3(csa, query, error)

  len = query.length

  l1 = len/3
  e1 = error/3
  q1 = query[0..l1-1]

  l2 = (len-l1)/2
  e2 = (error-e1)/2
  q2 = query[l1..l1+l2-1]

  l3 = len - l1 - l1
  e3 = error - e1 - e2
  q3 = query[l1+l2..len-1]

  l12 = l1+l2
  q12 = query[0..l1+l2-1]

  l23 = l2+l3
  q23 = query[l1..len-1]



  print("l1 = #{l1} e1 = #{e1} l2 = #{l2} e2 = #{e2} l3 = #{l3} e3 = #{e3}\n")

  list = []
  range = [0,csa.length()]

  d = Array.new(l1+1){|i| i}
  list1 = approx_sub(csa, q1, e1, d, range, 0, 0, 0)
  d1 = Array.new(l23+1)
  list1.each {|e, dofs, range2|
    for i in 0..l23 do d1[i] = e+i end
    tmp = approx_sub2(csa, q23, error, d1, range2, 0, dofs, 1)
    if tmp != [] then list.concat(tmp) end
  }

  d = Array.new(l3+1){|i| i}
  list1 = approx_sub(csa, q3, e3, d, range, 0, 0, 0)
  d1 = Array.new(l12+1)
  list1.each {|e, dofs, range2|
    for i in 0..l12 do d1[i] = e+i end
    tmp = approx_sub2(csa, q1, error, d1, range2, 0, dofs, 0)
    if tmp != [] then list.concat(tmp) end
  }

  list2 = []
  d = Array.new(l2+1){|i| i}
  list1 = approx_sub2(csa, q2, e2, d, range, 0, 0, 0)
  d1 = Array.new(l1+1)
  list1.each {|e, dofs, range2|
    for i in 0..l1 do d1[i] = e+i end
    tmp = approx_sub2(csa, q1, e2*2, d1, range2, 0, dofs, 0)
    if tmp != [] then list2.concat(tmp) end
  }
  d1 = Array.new(l3+1)
  list2.each {|e, dofs, range2|
    for i in 0..l3 do d1[i] = e+i end
    tmp = approx_sub2(csa, q3, error, d1, range2, 0, dofs, 1)
    if tmp != [] then list.concat(tmp) end
  }


  ans = []
  list.each {|e, len, range2|
    ee = dp(query,csa.substring(range2[0], len),0)
    ans << [ee, len, range2]
  }
#  ans.sort!{|a,b| [a[0],a[2],a[1]] <=> [b[0],b[2],b[1]]}
#  ans.uniq!
  ans = merge(csa, error, ans)

  return ans
end


csa = CSA.new(ARGV[0],ARGV[1])

#n = csa.length()
#pos = 100000
pos = 200000

while 1 do
  print("input score len ")
  buf = STDIN.gets().chomp.split(" ")
  k = buf[0].to_i
  len = buf[1].to_i
  print("k = #{k} len = #{len}\n")
  query = csa.text(pos, pos+len-1)
  print("query = ", query, "\n")

=begin
  $nnodes = 0
  ans = nil
  puts Benchmark.measure{
    ans = approx(csa, query, k)
  }
  print("#visited nodes = #{$nnodes}\n")

=end

=begin
  ans.each {|e,len,x|
    print("[#{x[0]}, #{x[1]}] error=#{e} \n")
    dp(query,csa.substring(x[0],len),1)
#    for i in x[0]..x[1]
#      print("#{csa.lookup(i)} ")
#    end
#    puts
  }
=end

  $nnodes = 0
  ans2 = nil
  puts Benchmark.measure{
    ans2 = approx2(csa, query, k)
  }
  print("#visited nodes = #{$nnodes}\n")

begin
  ans2.each {|e,len,x|
    print("[#{x[0]}, #{x[1]}] error=#{e} \n")
    dp(query,csa.substring(x[0],len),1)
#    for i in x[0]..x[1]
#      print("#{csa.lookup(i)} ")
#    end
#    puts
  }
end

=begin
  $nnodes = 0
  ans3 = nil
  puts Benchmark.measure{
    ans3 = approx3(csa, query, k)
  }
  print("#visited nodes = #{$nnodes}\n")
=end

=begin
  ans3.each {|e,len,x|
    print("[#{x[0]}, #{x[1]}] error=#{e} \n")
    dp(query,csa.substring(x[0],len),1)
#    for i in x[0]..x[1]
#      print("#{csa.lookup(i)} ")
#    end
#    puts
  }
=end


=begin
  while ans2 != [] || ans3 != [] do
    if ans2 != [] then
      l1 = ans2[0]
      e1 = l1[0];  len1 = l1[1];  x1 = l1[2]
      print("ans2 ")
      print("[#{x1[0]}, #{x1[1]}] error=#{e1} \n")
      dp(query,csa.substring(x1[0],len1),1)
      ans2.shift
    end
    if ans3 != [] then
      l2 = ans3[0]
      e2 = l2[0];  len2 = l2[1];  x2 = l2[2]
      print("ans3 ")
      print("[#{x2[0]}, #{x2[1]}] error=#{e2} \n")
      dp(query,csa.substring(x2[0],len2),1)
      ans3.shift
    end
  end
=end

end
