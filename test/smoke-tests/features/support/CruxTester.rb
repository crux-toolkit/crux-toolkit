require "CSV"
require "fileutils"
require "open3"
require "set"

class String
  def nan?
    self !~ /^\s*[+-]?((\d+_?)*\d+(\.(\d+_?)*\d+)?|\.(\d+_?)*\d+)(\s*|([eE][+-]?(\d+_?)*\d+)\s*)$/
  end
end

class CruxTester
  def initialize(path)
    @crux_path = path
    @crux_args = Array.new
  end

  def set_test_name(name) @crux_test_name = name end
  def add_arg(arg) @crux_args.push(arg) end

  def exec(cmd)
    if @crux_test_name == nil
      raise("set_test_name must be called before exec")
    end
    unless File.executable?(@crux_path)
      raise(@crux_path + " cannot be executed")
    end
    Open3.popen3(@crux_path + " " + cmd + " " + @crux_args.join(" ")) do | stdin, stdout, stderr, thread |
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
    same = FileUtils.cmp(expected, actual)
    copyObserved(actual, expected, same, write_observed)
    return same
  end

  def cmpTableWithPrecision(expected, actual, tolerance)

    # Compare the fields of a tab delimited table. Allow floating point
    # values to differ, as long as the relative error is less than
    # the given tolerance

    tolerance = tolerance.to_f

    unless File.readable?(expected)
      raise("cannot read file '" + expected + "'")
    end
    unless File.readable?(actual)
      raise("cannot read file '" + actual + "'")
    end
    expected_text = CSV.open(expected, {:col_sep => "\t"}).read
    actual_text = CSV.open(actual, {:col_sep => "\t"}).read
    # Compare tables line by line
    for i in 0..expected_text.size
      if expected_text[i] == actual_text[i]
        same = true
      else
        # Lines differ, compare value by value
        for j in 0..expected_text[i].size
          if expected_text[i][j] == actual_text[i][j]
            same = true
          else
            # Values differ can we compare as floats?
            if expected_text[i][j].nan? or actual_text[i][j].nan?
              # Apparently they're not floats
              same = false
              break
            else
              # They're floats, do they match within the relative error?
              expected_value = expected_text[i][j].to_f
              actual_value = actual_text[i][j].to_f
              relative_difference = ((actual_value - expected_value) /  expected_value).abs
              if  relative_difference <= tolerance
                same = true
              else
                same = false
                break
              end
            end
          end
        end
        if not same then break end
      end
    end
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

