% Script to check for the missing ground motions based on the flatfile
clc
clear


a = readtable('Flatfile-NGA-West2-Reduced.csv');
RSN = a.record_sequence_number_NGA;

% get the tokens of the files in the directory
b = dir('NGA-West2-Reduced/RSN*_1.AT2');
c = dir('NGA-West2-Reduced/RSN*_2.AT2');
d = dir('NGA-West2-Reduced/RSN*_3.AT2');
names1={b.name};
names2={c.name};
names3={d.name};
names1 = strrep(names1,'_1.AT2','') ; % remove extension
names2 = strrep(names2,'_2.AT2','') ; % remove extension
names3 = strrep(names3,'_3.AT2','') ; % remove extension
ids1 = str2double(strrep(names1,'RSN','')) ; % remove extension
ids2 = str2double(strrep(names2,'RSN','')) ; % remove extension
ids3 = str2double(strrep(names3,'RSN','')) ; % remove extension

% find missing gms
miss1 = setdiff(RSN,ids1)
miss2 = setdiff(RSN,ids2)
miss3 = setdiff(RSN,ids3)