require "rspec/expectations"

Given /^the working directory is (.*)$/ do | dir |
  @dir = dir
  Dir.chdir(dir)
end

Given /^the path to [Cc]rux is (.*)$/ do | path |
  @tester = CruxTester.new(path)
end

Given /^I want to run a test named (.*)$/ do | test_name |
  @tester.set_test_name(test_name)
end

Given /^I pass the arguments? (.*)$/ do | arg |
  @tester.add_arg(arg)
end

When /^I run ([^\s]+)$/ do | cmd |
  @last_ret = @tester.exec(cmd)
  @tester.test_done()
end

When /^I run ([^\s]+) as an intermediate step$/ do | cmd |
  @last_ret = @tester.exec(cmd)
end

When /^I ignore lines matching the pattern: \/(.*)\/$/ do | pattern |
  @tester.add_ignore_pattern(pattern)
end

Then /^the return value should be (-?[0-9]+)$/ do | ret |
  expect(@last_ret).to eq(ret.to_i)
end

Then /^(.*) should match (.*)$/ do | actual, expected |
  expect(@tester.cmp(expected, actual)).to be true
end

Then /^(.*) should contain the same lines as (.*)$/ do | actual, expected |
  expect(@tester.cmpUnordered(expected, actual)).to be true
end

