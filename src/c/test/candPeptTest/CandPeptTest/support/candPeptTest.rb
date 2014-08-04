class CandPeptTest
	def initialize()
		@crux_args = Array.new;
	end
	
	def crux_path(n)
		@crux_path = n;
	end
	
	def comet_path(n)
		@comet_path = n;
	end

	def fasta(n)
		@fasta_list = n.split(",");
	end

	def msms(n)
		@msms_list = n.split(",");
	end
	
	def mass(n)
		@mass_list = n.split(",");
	end
	
	def mass_type(n)
		@mass_type_list = n.split(",");
	end
	
	def missed_cleavages(n)
		@missed_cleavages_list = n.split(",");
	end

	def add_args(n)
		@crux_args.push(n);
	end

	def runCommand()
		mtype = 0;
		for mass in @mass_list do
			for mass_type in @mass_type_list do
				for fasta in @fasta_list do
					for msms in @msms_list do
						for mc in @missed_cleavages_list do
							#print parameter file
							File.open("crux.prm", 'w') { |file|
								file.puts("# comet_version 2014.01 rev. 0");
								file.puts("database_name = "+fasta);
								file.puts("precursor-window="+mass);
								file.puts("peptide_mass_tolerance="+mass);
								file.puts("precursor-window-type="+mass_type);
								if (mass_type.eql? "ppm")
									file.puts("peptide_mass_units=2");
									mtype = 0;
								elsif (mass_type.eql? "mass")
									file.puts("peptide_mass_units=0");
									mtype = 1;
								end		
								file.puts("missed-cleavages="+mc);
								file.puts("allowed_missed_cleavage="+mc);
							}
							File.open("crux.prm", 'a') { |file| @crux_args.each{|x| file.puts(x) } }
							#execute tide index
							cmd = @crux_path + " tide-index " + fasta +" crux-output" +
											" --parameter-file crux.prm --overwrite T ";
							puts cmd;
							system(cmd);
							#execute tide search
							cmd = @crux_path + " tide-search --parameter-file crux.prm" + 
											 " --overwrite T " + msms +" crux-output";
							system(cmd);
							#execute comet
							cmd = @crux_path + " comet " + msms   +  " " + fasta +
												" --parameter-file crux.prm --overwrite T";
							#cmd = @comet_path + " -Pcrux.prm " + msms;						
							system(cmd);
							return false if compareCandidatePeptideSets(mass, mtype) == false;
						end
					end
				end
			end
		end
		return true;
	end

	def compareCandidatePeptideSets(mass, type)
		mass = mass.to_f;
		@Table1 = Hash.new();
		#Open and parse the first (TIDE-SEARCH) output file
		file = File.new("crux-output/tide-search.target.txt", "r");
		line = file.gets;  # drop the first header line
		while (row = file.gets)
			record = row.split("\t");
			next if (record.length() < 10);
			@PSM = record[0]+record[1]+record[9]; 
			@Table1[@PSM] = 1;
			if (type == 0 && ( ((record[4].to_f-record[3].to_f)/record[4].to_f).abs > mass*0.9 ) ) #mass tolerance type is ppm
				@Table1[@PSM] = 2;
			end
			if (type == 1 && ( (record[4].to_f-record[3].to_f).abs > mass*0.9 ) ) #mass tolerance type is 'mass'
				@Table1[@PSM] = 2;
			end			
		end
		file.close;

		#open and parse the second (COMET) file
		file = File.new("crux-output/comet.target.txt", "r");
		line = file.gets;  # drop the first header line
#		line = file.gets;  # drop the first header line
		while (row = file.gets)
			record = row.split("\t"); 
			next if (record.length() < 10);
			@PSM = record[0]+record[1]+record[13]; 
			if (@Table1[@PSM] == nil)
				if (type == 0 && ( ((record[3].to_f - record[2].to_f)/record[3].to_f).abs > mass*0.9 )  )  #mass tolerance type is ppm
					return false;
				end
				if (type == 1 && ( (record[4].to_f - record[3].to_f).abs > mass*0.9 ) ) #mass tolerance type is 'mass'
					return false;
				end
			else
				@Table1[@PSM] += 1;
			end
		end
		file.close;
		#find any element in Table1 that has not been accessed.
		if (@Table1.has_value?(1) == true)
			return false;
		end
		return true;
	end
end
