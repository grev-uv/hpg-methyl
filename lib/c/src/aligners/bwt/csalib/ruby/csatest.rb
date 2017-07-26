require 'csa'

csa = CSA.new('E.coli.idx','E.coli.bwd')
for i in 0..9
  j = csa.lookup(i)
  print("SA[", i, "] = ", j, "\n")
end

while 1 do
  print("input score len ")
  key = gets().chomp
  puts "key => #{key.inspect}"
  break if key.empty?

  range = csa.search(key)
  if range != nil then
    puts "[#{range[0]}, #{range[1]}]"
    for i in range[0]..(range[0]+20)
      j = csa.lookup(i)
      puts "#{j}: #{csa.text(j-20, j+20)}"
      puts "#{j}: #{csa.substring(i, 20)}"
    end
  else
    puts "not found"
  end
end
