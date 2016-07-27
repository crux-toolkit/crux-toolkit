require "fileutils"
require "open3"
require "set"

class CruxTester
  def initialize(path)
    @crux_path = path
    @crux_args = Array.new
    @ignore_patterns = Array.new
  end

  def set_test_name(name) @crux_test_name = name end
  def add_arg(arg) @crux_args.push(arg) end
  def add_ignore_pattern(pattern) @ignore_patterns.push(Regexp.new(pattern)) end

  def exec(cmd)
    if @crux_test_name == nil
      raise("set_test_name must be called before exec")
    end
    unless File.executable?(@crux_path)
      raise(@crux_path + " cannot be executed")
    end
    Open3.popen3(@crux_path + " " + cmd + " --no-analytics T " + @crux_args.join(" ")) do | stdin, stdout, stderr, thread |
      @last_stdout = ""
      while line = stdout.gets
        @last_stdout << line
      end
      if @last_stdout.empty?
        @last_stdout = nil
      end
      while line = stderr.gets
        #puts("2>" + line)
      end
      @crux_args = Array.new
      return 0
      #TODO return thread.value.exitstatus
    end
  end

  def test_done()
    @crux_test_name = nil
  end

  def cmp(expected, actual, write_observed = 1)
    if actual == "stdout"
      same = @last_stdout == File.read(expected)
      writeObserved(@last_stdout, expected, same, write_observed)
      return same
    end

    if File.directory?(expected) && File.directory?(actual)
      Dir.foreach(expected) do | dirent |
        next if dirent == "." or dirent == ".."
        unless cmp(expected + "/" + dirent, actual + "/" + dirent)
          return false
        end
      end
      return true
    end

    unless File.readable?(expected)
      raise("cannot read file '" + expected + "'")
    end
    unless File.readable?(actual)
      raise("cannot read file '" + actual + "'")
    end

    same = true
    file_expected = File.open(expected)
    file_actual = File.open(actual)
    file_expected.each.zip(file_actual.each).each do | line_expected, line_actual |
      if not line_expected == line_actual
        ignore = false
        @ignore_patterns.each do | pattern |
          if pattern.match(line_expected) != nil && pattern.match(line_actual) != nil
            ignore = true
            break
          end
        end
        if not ignore
          same = false
          break
        end
      end
    end
    file_expected.close
    file_actual.close
    copyObserved(actual, expected, same, write_observed)
    return same
  end

  # Compare files where line order does not matter
  def cmpUnordered(expected, actual, write_observed = 1)
    if File.directory?(expected) && File.directory?(actual)
      Dir.foreach(expected) do | dirent |
        next if dirent == "." or dirent == ".."
        unless cmpUnordered(expected + "/" + dirent, actual + "/" + dirent)
          return false
        end
      end
      return true
    end

    unless File.readable?(expected)
      raise("cannot read file '" + expected + "'")
    end
    unless File.readable?(actual)
      raise("cannot read file '" + actual + "'")
    end
    expected_lines = 0
    actual_lines = 0
    only_expected = Set.new
    only_actual = Set.new
    file_expected = File.open(expected)
    file_actual = File.open(actual)
    file_expected.each.zip(file_actual.each).each do | line_expected, line_actual |
      if not line_expected == nil
        expected_lines += 1
      end
      if not line_actual == nil
        actual_lines += 1
      end
      if not expected_lines == actual_lines
        break
      elsif not line_expected == line_actual
        if only_expected.include?(line_actual)
          only_expected.delete(line_actual)
        else
          only_actual.add(line_actual)
        end
        if only_actual.include?(line_expected)
          only_actual.delete(line_expected)
        else
          only_expected.add(line_expected)
        end
      end
    end
    file_expected.close
    file_actual.close
    same = only_expected.empty?() && only_actual.empty?() && expected_lines == actual_lines
    copyObserved(actual, expected, same, write_observed)
    return same
  end

  def writeObserved(actual_content, expected_filename, same, write_observed)
    if write_observed == 1
      observed = expected_filename + ".observed"
      if not same
        File.open(observed, "w") { |file| file.write(actual_content) }
      elsif File.file?(observed)
        FileUtils.rm(observed)
      end
    end
  end

  def copyObserved(actual_filename, expected_filename, same, write_observed)
    if write_observed == 1
      observed = expected_filename + ".observed"
      if not same
        FileUtils.cp(actual_filename, expected_filename + ".observed")
      elsif File.file?(observed)
        FileUtils.rm(observed)
      end
    end
  end
end

