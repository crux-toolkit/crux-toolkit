require "fileutils"
require "open3"

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
      @crux_test_name = nil
      @crux_args = Array.new
      return 0
      #TODO return thread.value.exitstatus
    end
  end

  def self.cmp(expected, actual)
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
    return FileUtils.cmp(expected, actual)
  end
end

