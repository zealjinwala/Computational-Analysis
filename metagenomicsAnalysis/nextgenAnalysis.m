%% Metagenomics Analysis
% By Zeal Jinwala

DATADIR = [bmes.datadir '/hwmetagenome/']; if ~bmes.isfolder(DATADIR); mkdir(DATADIR); end
BWAEXE = bmes.bwaexe; 

%% Download NGS data
fastqfile = 'http://sacan.biomed.drexel.edu/ftp/SRR3656745_pass.randsample.select.102.fastq';
if ~isempty(regexp(fastqfile,'^(https?://|ftps?://)','once'))
	fastqfile = bmes.downloadurl(fastqfile, DATADIR);
end
reads = fastqread(fastqfile);

%% Download the Genomes
% Olsenella uli. NC_014363.
OlsenallaUli = bmes.downloadurl('https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=302334808&extrafeat=null&conwithfeat=on',[DATADIR '/OlsenellaUli.NC_014363.1.fasta']);
% Segniliparus rotundu. NC_014168. Find the genome of this organism on NCBI and download it as a fasta file
Segniliparus = bmes.downloadurl('https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=296392440&',[DATADIR '/Segniliparus.NC_014168.1.fasta']);
% Escherichia coli K-12. NC_000913. Find the genome of this organism on NCBI and download it as a fasta file.
Escherichia = bmes.downloadurl('https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=556503834&',[DATADIR '/Escherichia.NC_000913.1.fasta']);
genomes = {OlsenallaUli Segniliparus Escherichia};

%% Index the genome files and Align sequence reads to the reference genomes
    % fastq file -> sequence reads files
    % fasta file -> genome files
    % using BWA to index genomes and aligning genomes to sequence reads

N = numel(genomes);
for i=1:N
    [~, name] = fileparts(genomes{i});
    % samfile name for each genome
    samfile = [fastqfile '.' name '.sam'];
    % indexing the genome 
    if ~bmes.isfileandnotempty([genomes{i} '.bwt'])
        cmd=[BWAEXE ' index "' genomes{i} '"' ];
        fprintf('Executing command: %s\n', cmd);
        system(cmd);
    end 
    
    % mapping the reads (fastq file) to the genome
    if ~bmes.isfileandnotempty(samfile)
         cmd=[BWAEXE ' mem "' genomes{i} '"' ' "' fastqfile '" > "' samfile '"'];
         fprintf('Executing command: %s\n', cmd);
         system(cmd);
    end
    
    % storing mapped data for all genomes
    sams{i} = samread(samfile);
end

%% Check the result
    % Check the result of the BWA to assign each short read to one of the genomes, or assign it "Unknown" if it cannot be mapped to any of the genomes.
    % If a short read does not map to any organism, increase the count for the Unknown group by 1. Not-aligned short reads have a ReferenceName '*' and a Position 0.
    % If a short read maps to only one organism, increase the count for that organism by 1.
    % If a short read maps to more than one organism, use a fractional increase in the organism counts. E.g., if a read is mapped to the first and third organisms only, increase the counts for the first and third organisms by 0.5, and do not increase the count for the second organism. E.g., if a read is mapped to all three organisms; increase each of their counts by (1/3).
    % Note: bwa may return more than one hit within a genome for a read. You need to make sure such hits are not counted multiple times for a genome.

% number of reads [rows] x number of genomes (including unknown) [cols]
mapMatrix = zeros(numel(reads), N+1); 
% store all reads from the fastq file in 1 cell array
readsSeqs = cellstr(char(reads.Sequence)); 

% For each genome
for i = 1:N
    % Sequence information for that genome
    genomeSeqs = cellstr(char(sams{i}.Sequence)); 
    % Storing genome lenghts to for that genome 
    gLengths(i) = numel(horzcat(sams{i}.Sequence)); 
    % for each sequence read in the fastqfile
    for j = 1:numel(readsSeqs)
			if j==4&&i==2
				2+2;
			end
        % If a short read maps to only one organism, increase the count for that organism by 1
         if sum(strcmpi(readsSeqs(j), genomeSeqs)) ~= 0
             mapMatrix(j,i) = 1;
         else
        % if a short read does not map to any organism, increase the count for the Unknown group by 1.
             mapMatrix(j,4) = mapMatrix(j,4) + 1;
         end 
    end 
end
return
% Checking bwa results and normalizing count: 
% If a short read maps to more than one organism, use a fractional increase in the organism counts. 
mapMatrix(mapMatrix(:,4) < 3 , 4) = 0;
mapMatrix(mapMatrix(:,4) == 3 , 4) = 1;
locs = find(sum(mapMatrix(:,1:3),2) > 1);
zeroCols = mapMatrix(:,1:3) == 0;
mapMatrix(locs,1:3) = repmat(1./(sum(mapMatrix(locs,1:3),2)),1,3);
mapMatrix(zeroCols) = 0;

%% Results
    % Show a bar graph of total read counts assigned to each organism and to the Unknown group. 
    % The x axis of the bar graph should show the organisms (including the Unknown group), 
    % and the y axis should show the count of reads mapped to each organism.

% get total number of reads mapped to each genome using normalized counts
totalCounts = sum(mapMatrix, 1);
bar(totalCounts)
hold on
legend = {'Olsenella uli','Segniliparus rotundu','Escherichia coli K-12', 'Unknown'};
set(gca,'Xtick',[1:4],'XtickLabel', legend, 'FontSize', 12);
xtickangle(30)
title('Number of Reads Mapped per Organism')
xlabel('Organisms')
ylabel('Number of Reads')
hold off

    % Show the percent abundances of the species (excluding the Unknown group) 
    % as a bar graph. Label the bars with the species names. 
    % Make sure you normalize the counts of the reads assigned to each organism by the genome size of that organism. 
    % The percentages of the three organisms should add up to 100%.
    % normalize counts by the length by diving total count by the length of each genome
normVals = totalCounts(1:3) ./ gLengths;
% calculate percent abundances
perc = (normVals ./ sum(normVals)) * 100;
bar(perc)
hold on
set(gca,'Xtick',[1:3], 'XtickLabel', legend(1:3), 'FontSize', 12);
xtickangle(30)
title('Percent abundances per Organism')
xlabel('Organisms')
ylabel('Percent abundances')
hold off