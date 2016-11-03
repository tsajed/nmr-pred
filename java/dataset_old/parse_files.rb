Dir.foreach('.') do |item|
  next if item == '.' or item == '..'
  if /21\-53/.match(item)
    name = item.split(' ')[0]
    File.open(name, 'w') do |f|
      File.open(item, 'r') do |i|
        f.write(i.read)
      end
    end
  end
end