use strict;
use warnings;
use File::Basename;
use XML::Simple;
use Data::Dumper;


my @dirs=glob("*");                                                      #获取所有目录
open(WF,">clinicalpractice.xls") or die $!;
print WF "Id\tfutime\tfustat\tgender\tage\tstage\tT\tM\tN\tChild_pugh_classification\tGrade\tRace\tAFP_at_procurement\n";
foreach my $dir(@dirs){                                                  #以下为爬虫处理
	if(-d $dir){                                                         #如果为目录下文件则处理
	  opendir(RD,"$dir") or die $!;                                      #打开文件
	  while(my $xmlfile=readdir(RD)){
	  	if($xmlfile=~/\.xml$/){                                          #判断是否xml文件，是则读取
	  		#print "$dir\\$xmlfile\n";
				my $userxs = XML::Simple->new(KeyAttr => "name");
				my $userxml = $userxs->XMLin("$dir\\$xmlfile");
				my $disease_code=$userxml->{'admin:admin'}{'admin:disease_code'}{'content'};   #提取文件信息，get disease code
				my $disease_code_lc=lc($disease_code);
				my $patient_key=$disease_code_lc . ':patient';                                #ucec:patient
				my $follow_key=$disease_code_lc . ':follow_ups';
				my $patient_barcode=$userxml->{$patient_key}{'shared:bcr_patient_barcode'}{'content'};  #TCGA-AX-A1CJ
				my $gender=$userxml->{$patient_key}{'shared:gender'}{'content'};#性别
				my $age=$userxml->{$patient_key}{'clin_shared:age_at_initial_pathologic_diagnosis'}{'content'};#年龄
				my $race=$userxml->{$patient_key}{'clin_shared:race_list'}{'clin_shared:race'}{'content'};#人种
				my $grade=$userxml->{$patient_key}{'shared:neoplasm_histologic_grade'}{'content'};  #G1/G2/G3				
			    my $pathologic_stage=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:pathologic_stage'}{'content'};  #stage I
				my $pathologic_T=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:pathologic_categories'}{'shared_stage:pathologic_T'}{'content'};
				my $pathologic_M=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:pathologic_categories'}{'shared_stage:pathologic_M'}{'content'};
				my $pathologic_N=$userxml->{$patient_key}{'shared_stage:stage_event'}{'shared_stage:tnm_categories'}{'shared_stage:pathologic_categories'}{'shared_stage:pathologic_N'}{'content'};
				my $child_pugh_classification=$userxml->{$patient_key}{'chol_lihc_shared:child_pugh_classification_grade'}{'content'};
				my $AFP=$userxml->{$patient_key}{'chol_lihc_shared:fetoprotein_outcome_value'}{'content'};
				
				
				$gender=(defined $gender)?$gender:"unknow";
				$age=(defined $age)?$age:"unknow";
				$race=(defined $race)?$race:"unknow";
				$grade=(defined $grade)?$grade:"unknow";
				$pathologic_stage=(defined $pathologic_stage)?$pathologic_stage:"unknow";
				$pathologic_T=(defined $pathologic_T)?$pathologic_T:"unknow";
				$pathologic_M=(defined $pathologic_M)?$pathologic_M:"unknow";
				$pathologic_N=(defined $pathologic_N)?$pathologic_N:"unknow";
				$child_pugh_classification=(defined $child_pugh_classification)?$child_pugh_classification:"unknow";
				$AFP=(defined $AFP)?$AFP:"unknow";
				
								
				my $survivalTime="";
				my $vital_status=$userxml->{$patient_key}{'clin_shared:vital_status'}{'content'};
				my $followup=$userxml->{$patient_key}{'clin_shared:days_to_last_followup'}{'content'};
				my $death=$userxml->{$patient_key}{'clin_shared:days_to_death'}{'content'};
				if($vital_status eq 'Alive'){
					$survivalTime="$followup\t0";
				}
				else{
					$survivalTime="$death\t1";
				}
				
				for my $i(keys %{$userxml->{$patient_key}{$follow_key}}){
					my @survivalArr=split(/\t/,$survivalTime);
					eval{
						$followup=$userxml->{$patient_key}{$follow_key}{$i}{'clin_shared:days_to_last_followup'}{'content'};
						$vital_status=$userxml->{$patient_key}{$follow_key}{$i}{'clin_shared:vital_status'}{'content'};
						$death=$userxml->{$patient_key}{$follow_key}{$i}{'clin_shared:days_to_death'}{'content'};
				  };
				  if($@){
					  $followup=$userxml->{$patient_key}{$follow_key}{$i}[0]{'clin_shared:days_to_last_followup'}{'content'};
						$vital_status=$userxml->{$patient_key}{$follow_key}{$i}[0]{'clin_shared:vital_status'}{'content'};
						$death=$userxml->{$patient_key}{$follow_key}{$i}[0]{'clin_shared:days_to_death'}{'content'};
				  }
					if($vital_status eq 'Alive'){
						if($followup>$survivalArr[0]){
					    $survivalTime="$followup\t0";
					  }
				  }
				  else{
				  	if($death>$survivalArr[0]){
					    $survivalTime="$death\t1";
					  }
				  }
				}
				print WF "$patient_barcode\t$survivalTime\t$gender\t$age\t$pathologic_stage\t$pathologic_T\t$pathologic_M\t$pathologic_N\t$child_pugh_classification\t$grade\t$race\t$AFP\n";
			}                        #上一行将所有信息打印出来
		}
		close(RD);
	}
}
close(WF);

