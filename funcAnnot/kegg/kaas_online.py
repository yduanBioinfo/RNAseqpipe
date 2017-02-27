#!/usr/bin/env python
#coding=utf-8

'''
it's a script for convenient using kaas web tool

_author : DY
_version : 0.0.2
_python : 2.6

'''

import urllib2 as urr
import urllib, sys, re
from bs4 import BeautifulSoup
from poster.encode import multipart_encode
from poster.streaminghttp import register_openers
	
url = 'http://www.genome.jp/kaas-bin/kaas_main'

#parameters
prog = 'BLAST' #BLAST/GHOSTX/GHOSTZ (amino acid query only)
peptide2 = 'n' #'n'/''
mail = 'username@example.com'
#org_list = 'hsa, dme, cel, ath, sce, cho, eco, nme, hpy, rpr, bsu, lla, cac, mge, mtu, ctr, bbu, syn, bth, dra, aae, mja, ape, ptr'
org_list = 'hsa, dre'
way = 'b'#b/s

#needn't alter
mycontinue = 1
uptype = 'q_text'
text = ''
file = ''
qname = 'query_name'
mode = 'compute'
org_ids = ['hsa', 'ptr', 'pps', 'ggo', 'pon', 'nle', 'mcc', 'mcf', 'rro', 'cjc', 'mmu', 'rno', 'cge', 'ngi', 'hgl', 'ocu', 'tup', 'cfa', 'aml', 'umr', 'fca', 'ptg', 'bta', 'bom', 'phd', 'chx', 'oas', 'ssc', 'cfr', 'bacu', 'lve', 'ecb', 'myb', 'myd', 'pale', 'mdo', 'shr', 'oaa', 'gga', 'mgp', 'apla', 'tgu', 'gfr', 'fab', 'phi', 'ccw', 'fpg', 'fch', 'clv', 'asn', 'amj', 'pss', 'cmy', 'acs', 'pbi', 'xla', 'xtr', 'dre', 'tru', 'mze', 'ola', 'xma', 'lcm', 'cmk', 'bfo', 'cin', 'spu', 'dme', 'dpo', 'dan', 'aga', 'aag', 'cqu', 'ame', 'nvi', 'tca', 'bmor', 'api', 'isc', 'cel', 'cbr', 'bmy', 'loa', 'tsp', 'hro', 'lgi', 'crg', 'smm', 'nve', 'hmg', 'tad', 'aqu', 'ath', 'aly', 'crb', 'brp', 'cit', 'tcc', 'gmx', 'fve', 'csv', 'vvi', 'sly', 'osa', 'olu', 'ota', 'mis', 'cme', 'gsl', 'sce', 'ago', 'kla', 'vpo', 'zro', 'cgr', 'ncs', 'tpf', 'ppa', 'dha', 'pic', 'pgu', 'lel', 'cal', 'yli', 'clu', 'ncr', 'mgr', 'fgr', 'nhe', 'maw', 'ssl', 'bfu', 'ani', 'afm', 'aor', 'ang', 'nfi', 'pcs', 'cim', 'cpw', 'pbl', 'ure', 'pno', 'bze', 'tml', 'spo', 'cne', 'cgi', 'ppl', 'mpr', 'scm', 'uma', 'mgl', 'ecu', 'nce', 'mbr', 'ddi', 'dfa', 'ehi', 'edi', 'acan', 'pfa', 'pyo', 'pkn', 'tan', 'tpv', 'bbo', 'cpv', 'cho', 'tgo', 'tet', 'ptm', 'pti', 'tps', 'ehx', 'gtt', 'tbr', 'tcr', 'lma', 'lif', 'lbz', 'ngr', 'tva', 'gla', 'eco', 'ecj', 'ecd', 'ebw', 'ece', 'ecs', 'ecf', 'etw', 'eoj', 'eoi', 'eoh', 'ecg', 'eok', 'ecc', 'ecp', 'eci', 'ecv', 'ecx', 'ecw', 'ecm', 'ecy', 'ecr', 'ecq', 'eck', 'ect', 'eum', 'ecz', 'esm', 'ecl', 'ebr', 'ebd', 'efe', 'sty', 'stt', 'sent', 'stm', 'setu', 'spt', 'sek', 'spq', 'sei', 'sec', 'seh', 'see', 'sew', 'sea', 'sed', 'seg', 'set', 'senj', 'ses', 'sbg', 'ype', 'ypk', 'ypa', 'ypn', 'ypm', 'ypp', 'ypg', 'ypz', 'yps', 'ypi', 'ypy', 'ypb', 'yen', 'yep', 'sfl', 'sfx', 'sfv', 'ssn', 'sbo', 'sbc', 'sdy', 'eca', 'pct', 'pwa', 'eta', 'epy', 'eam', 'eay', 'ebi', 'plu', 'pay', 'buc', 'bap', 'bau', 'bas', 'bab', 'bcc', 'wbr', 'sgl', 'ent', 'enc', 'esc', 'esa', 'ctu', 'kpn', 'kpu', 'kpm', 'kpe', 'kva', 'kox', 'cko', 'cro', 'spe', 'pmr', 'eic', 'etr', 'bfl', 'bpn', 'hde', 'dda', 'dze', 'ddc', 'xbo', 'pam', 'rip', 'rah', 'men', 'hin', 'hit', 'hip', 'hiq', 'hif', 'hiu', 'hdu', 'hap', 'hpr', 'hso', 'hsm', 'pmu', 'msu', 'apl', 'apj', 'apa', 'asu', 'aap', 'aat', 'gan', 'xfa', 'xft', 'xfm', 'xfn', 'xcc', 'xcb', 'xca', 'xcv', 'xac', 'xoo', 'xom', 'xop', 'xal', 'sml', 'smt', 'psu', 'fau', 'vch', 'vcj', 'vco', 'vcm', 'vvu', 'vvy', 'vvm', 'vpa', 'vha', 'vsp', 'vex', 'van', 'vfi', 'vfm', 'vsa', 'ppr', 'pae', 'pau', 'pap', 'pag', 'ppu', 'ppf', 'ppg', 'ppw', 'pst', 'psb', 'psp', 'pfl', 'pfo', 'pfs', 'pen', 'pmy', 'psa', 'avn', 'par', 'pcr', 'prw', 'acb', 'abm', 'aby', 'abc', 'abn', 'abb', 'acd', 'aci', 'mct', 'son', 'sdn', 'sfr', 'saz', 'sbl', 'sbm', 'sbn', 'sbp', 'slo', 'spc', 'sse', 'spl', 'she', 'shm', 'shn', 'shw', 'shl', 'swd', 'swp', 'svo', 'ilo', 'cps', 'pha', 'pat', 'maq', 'amc', 'gag', 'pin', 'fbl', 'cja', 'sde', 'ttu', 'cbu', 'cbs', 'cbd', 'cbg', 'cbc', 'lpn', 'lpf', 'lpp', 'lpc', 'lpa', 'llo', 'mca', 'mmt', 'mah', 'ftu', 'ftf', 'ftw', 'ftl', 'fth', 'fta', 'ftm', 'ftn', 'fph', 'tcx', 'tcy', 'mej', 'noc', 'nhl', 'nwa', 'alv', 'aeh', 'hha', 'tgr', 'tkm', 'hna', 'hch', 'csa', 'hel', 'crp', 'abo', 'kko', 'mmw', 'aha', 'asa', 'tau', 'oce', 'dno', 'bci', 'rma', 'vok', 'gpb', 'nma', 'nme', 'nmc', 'nmn', 'nmi', 'ngo', 'ngk', 'cvi', 'lhk', 'pse', 'rso', 'rsc', 'rsl', 'rpi', 'rpf', 'reu', 'reh', 'rme', 'cti', 'bma', 'bmv', 'bml', 'bmn', 'bps', 'bpm', 'bpl', 'bpd', 'bpr', 'bte', 'bvi', 'bur', 'bcn', 'bch', 'bcm', 'bcj', 'bam', 'bac', 'bmu', 'bmj', 'bxe', 'bph', 'bpy', 'bgl', 'bug', 'bge', 'bgf', 'brh', 'pnu', 'pne', 'bpe', 'bpa', 'bbr', 'bpt', 'bav', 'axy', 'teq', 'put', 'aka', 'rfr', 'pol', 'pna', 'aav', 'ajs', 'dia', 'aaa', 'vei', 'dac', 'del', 'vap', 'vpe', 'ctt', 'adn', 'adk', 'rta', 'mpt', 'har', 'mms', 'hse', 'zin', 'cfu', 'lch', 'tin', 'neu', 'net', 'nmu', 'eba', 'azo', 'dar', 'tmz', 'dsu', 'tbd', 'mfa', 'mmb', 'meh', 'mei', 'app', 'tpn', 'slt', 'gca', 'hpy', 'hpj', 'hpa', 'hps', 'hpg', 'hpp', 'hpb', 'hpl', 'hhe', 'hac', 'hms', 'wsu', 'tdn', 'sku', 'cje', 'cjj', 'cju', 'cjr', 'cjd', 'cff', 'ccv', 'cha', 'cco', 'cla', 'abu', 'ant', 'sdl', 'nsa', 'nis', 'sun', 'nam', 'gsu', 'gme', 'gur', 'glo', 'gbm', 'geo', 'gem', 'pca', 'ppd', 'dvu', 'dvl', 'dvm', 'dde', 'dds', 'dma', 'dsa', 'lip', 'dba', 'drt', 'bba', 'bmx', 'dps', 'dak', 'dpr', 'dol', 'dal', 'dat', 'ade', 'acp', 'afw', 'ank', 'mxa', 'mfu', 'sur', 'scl', 'hoh', 'sat', 'dao', 'sfu', 'dbr', 'hmr', 'rpr', 'rty', 'rcm', 'rbe', 'rbo', 'rco', 'rfe', 'rak', 'rri', 'rrj', 'rms', 'rpk', 'raf', 'ots', 'ott', 'wol', 'wri', 'wpi', 'wbm', 'ama', 'amf', 'acn', 'aph', 'eru', 'erw', 'erg', 'ecn', 'ech', 'nse', 'nri', 'mmn', 'mlo', 'mes', 'pla', 'sme', 'smd', 'rhi', 'atu', 'ara', 'avi', 'ret', 'rec', 'rle', 'rlt', 'rlg', 'las', 'bme', 'bmi', 'bmf', 'bmb', 'bmc', 'baa', 'bms', 'bmt', 'bov', 'bcs', 'bmr', 'oan', 'bja', 'bra', 'bbt', 'rpa', 'rpb', 'rpc', 'rpd', 'rpe', 'rpt', 'nwi', 'nha', 'oca', 'bhe', 'bqu', 'bbk', 'btr', 'bgr', 'xau', 'azc', 'sno', 'mex', 'mea', 'mdi', 'mch', 'mrd', 'met', 'mpo', 'mno', 'bid', 'msl', 'hdn', 'rva', 'phl', 'hci', 'ccr', 'ccs', 'cak', 'cse', 'pzu', 'bsb', 'aex', 'sil', 'sit', 'rsp', 'rsh', 'rsq', 'rsk', 'rcp', 'jan', 'rde', 'pde', 'dsh', 'kvu', 'mmr', 'hne', 'hba', 'zmo', 'zmn', 'nar', 'sal', 'swi', 'sjp', 'eli', 'gox', 'gbe', 'acr', 'gdi', 'gdj', 'apt', 'rru', 'rce', 'mag', 'azl', 'pbr', 'mgm', 'pgv', 'pub', 'apb', 'afe', 'afr', 'bsu', 'bli', 'bld', 'bay', 'bao', 'bha', 'ban', 'bar', 'bat', 'bah', 'bai', 'bce', 'bca', 'bcz', 'bcr', 'bcb', 'bcu', 'bcg', 'bcq', 'bcx', 'bcy', 'btk', 'btl', 'btb', 'bwe', 'bcl', 'bpu', 'bpf', 'bmq', 'bmd', 'bck', 'oih', 'gka', 'gtn', 'gth', 'gwc', 'gyc', 'gct', 'afl', 'lsp', 'bse', 'sau', 'sav', 'saw', 'sah', 'saj', 'sam', 'sas', 'sar', 'sac', 'sax', 'saa', 'sao', 'sae', 'sad', 'sab', 'sep', 'ser', 'sha', 'ssp', 'sca', 'slg', 'ssd', 'mcl', 'lmo', 'lmf', 'lmh', 'lmc', 'lmn', 'lmy', 'lin', 'lwe', 'lsg', 'esi', 'eat', 'bbe', 'pjd', 'gym', 'aac', 'bts', 'lla', 'llk', 'llc', 'llm', 'spy', 'spz', 'spm', 'spg', 'sps', 'sph', 'spi', 'spj', 'spk', 'spf', 'spa', 'spb', 'soz', 'spn', 'spd', 'spr', 'spw', 'spx', 'sne', 'spv', 'snm', 'sjj', 'spp', 'snt', 'snc', 'sag', 'san', 'sak', 'smu', 'smc', 'stc', 'stl', 'ste', 'ssa', 'ssb', 'ssu', 'ssv', 'ssi', 'sss', 'sgo', 'seq', 'sez', 'seu', 'sub', 'sds', 'sga', 'sgg', 'smb', 'std', 'lpl', 'lpj', 'ljo', 'ljf', 'lac', 'lsa', 'lsl', 'ldb', 'lbu', 'lbr', 'lca', 'lcb', 'lga', 'lre', 'lrf', 'lhe', 'lfe', 'lrh', 'lrl', 'lcr', 'ppe', 'pce', 'efa', 'efc', 'mps', 'thl', 'ooe', 'lme', 'lci', 'lki', 'lec', 'lgs', 'wko', 'aur', 'crn', 'cac', 'cpe', 'cpf', 'cpr', 'ctc', 'cno', 'cbo', 'cba', 'cbh', 'cby', 'cbl', 'cbk', 'cbb', 'cbi', 'cbt', 'cbf', 'cbe', 'ckl', 'ckr', 'amt', 'aoe', 'asf', 'cth', 'cce', 'eha', 'ral', 'clo', 'bpb', 'cle', 'rho', 'cpy', 'ere', 'cdf', 'cdc', 'cdl', 'sth', 'swo', 'slp', 'dsy', 'dhd', 'drm', 'dae', 'pth', 'dau', 'tjr', 'sgy', 'hmo', 'eel', 'ova', 'tmr', 'say', 'tte', 'tex', 'tpd', 'tit', 'tmt', 'chy', 'tep', 'mta', 'adg', 'csc', 'ate', 'toc', 'ttm', 'cpo', 'tnr', 'mas', 'nth', 'hor', 'has', 'aar', 'fma', 'apr', 'vpr', 'ssg', 'med', 'afn', 'erh', 'mge', 'mpn', 'mpu', 'mpe', 'mga', 'mmy', 'mcp', 'mmo', 'mhy', 'mhj', 'mhp', 'msy', 'maa', 'mal', 'mat', 'mco', 'mho', 'mcd', 'uur', 'upa', 'uue', 'poy', 'ayw', 'pml', 'pal', 'acl', 'mfl', 'mtu', 'mtc', 'mra', 'mtf', 'mtb', 'mbo', 'mbb', 'mbt', 'mle', 'mlb', 'mpa', 'mav', 'msm', 'mul', 'mva', 'mgi', 'msp', 'mab', 'mmc', 'mkm', 'mjl', 'mmi', 'asd', 'cgl', 'cgb', 'cgt', 'cef', 'cdi', 'cjk', 'cur', 'car', 'ckp', 'nfa', 'rha', 'rer', 'rop', 'gbr', 'tpr', 'srt', 'sco', 'sma', 'sgr', 'scb', 'twh', 'tws', 'lxx', 'cmi', 'cms', 'mts', 'art', 'aau', 'ach', 'rsa', 'krh', 'mlu', 'rmu', 'bcv', 'bfa', 'jde', 'kse', 'xce', 'iva', 'ske', 'cfl', 'ica', 'pac', 'pak', 'pfr', 'mph', 'nca', 'kfl', 'tfu', 'nda', 'tcu', 'sro', 'fra', 'fre', 'fri', 'fal', 'ace', 'nml', 'gob', 'kra', 'sen', 'svi', 'amd', 'pdx', 'ami', 'stp', 'saq', 'mau', 'mil', 'vma', 'cai', 'sna', 'ahe', 'mcu', 'blo', 'blj', 'bln', 'blf', 'bll', 'blb', 'bad', 'bla', 'blc', 'blt', 'bde', 'gva', 'tbi', 'rxy', 'cwo', 'afo', 'ccu', 'shi', 'ele', 'apv', 'ols', 'cgo', 'ctr', 'cta', 'ctb', 'ctl', 'ctj', 'cmu', 'cpn', 'cpa', 'cpj', 'cpt', 'cca', 'cab', 'cfe', 'pcu', 'puv', 'wch', 'sng', 'ote', 'caa', 'amu', 'min', 'bbu', 'bbz', 'bga', 'baf', 'btu', 'bhr', 'bdu', 'bre', 'tpa', 'tpp', 'tde', 'ssm', 'lil', 'lic', 'lbj', 'lbl', 'lbi', 'lbf', 'bhy', 'brm', 'bpo', 'aba', 'aca', 'acm', 'tsa', 'sus', 'ctm', 'fsu', 'emi', 'rsd', 'fnu', 'ipo', 'lba', 'str', 'smf', 'gau', 'tai', 'aco', 'tli', 'rba', 'psl', 'plm', 'ipa', 'syn', 'syw', 'syc', 'syf', 'syd', 'sye', 'syg', 'syr', 'syx', 'syp', 'cya', 'cyb', 'tel', 'mar', 'cyt', 'cyp', 'cyc', 'cyn', 'cyh', 'cyj', 'amr', 'cyu', 'ter', 'gvi', 'ana', 'npu', 'ava', 'naz', 'pma', 'pmm', 'pmt', 'pmn', 'pmi', 'pmb', 'pmc', 'pmf', 'pmg', 'pmh', 'pmj', 'pme', 'bth', 'bfr', 'bfs', 'bvu', 'pgi', 'pgn', 'pdi', 'ppn', 'osp', 'aps', 'pru', 'sru', 'srm', 'rmr', 'cpi', 'phe', 'psn', 'shg', 'hhy', 'cmr', 'chu', 'dfe', 'sli', 'lby', 'rsi', 'mtt', 'aas', 'gfo', 'fjo', 'fps', 'fbr', 'coc', 'rbi', 'zpr', 'cat', 'ran', 'fbc', 'cao', 'cly', 'wvi', 'kdi', 'lan', 'mrs', 'fba', 'smg', 'sms', 'smh', 'sum', 'bbl', 'bpi', 'fte', 'cte', 'cpc', 'cch', 'cph', 'cpb', 'cli', 'pvi', 'plt', 'pph', 'paa', 'cts', 'det', 'deh', 'deb', 'dev', 'deg', 'dly', 'rrs', 'rca', 'cau', 'cag', 'chl', 'hau', 'tro', 'sti', 'atm', 'dra', 'dge', 'ddr', 'tra', 'tth', 'ttj', 'mrb', 'msv', 'opr', 'mhd', 'aae', 'hya', 'hth', 'tal', 'sul', 'saf', 'pmx', 'tam', 'dte', 'tma', 'tpt', 'trq', 'tna', 'tnp', 'tle', 'tme', 'taf', 'fno', 'pmo', 'kol', 'ccz', 'din', 'ddf', 'dap', 'cni', 'fsi', 'dth', 'dtu', 'tye', 'nde', 'tid', 'top', 'ttr', 'mja', 'mfe', 'mvu', 'mfs', 'mif', 'mig', 'mmp', 'mmq', 'mmx', 'mmz', 'mae', 'mvn', 'mvo', 'mok', 'mac', 'mba', 'mma', 'mbu', 'mmh', 'mev', 'mzh', 'mpy', 'mhz', 'mtp', 'mhu', 'mla', 'mem', 'mpi', 'mbn', 'mpl', 'mpd', 'rci', 'mth', 'mmg', 'mst', 'msi', 'mru', 'mel', 'mfv', 'mka', 'afu', 'apo', 'fpl', 'hal', 'hsl', 'hma', 'nph', 'hut', 'hmu', 'hje', 'hwa', 'hla', 'hvo', 'hbo', 'htu', 'nmg', 'hxa', 'nat', 'nge', 'hru', 'nou', 'sali', 'tac', 'tvo', 'pto', 'mer', 'pho', 'pab', 'pfu', 'tko', 'ton', 'tga', 'tsi', 'tba', 'abi', 'ape', 'smr', 'shc', 'iho', 'dka', 'dmu', 'tag', 'iag', 'thg', 'hbu', 'pfm', 'sso', 'sto', 'sai', 'sis', 'sia', 'sim', 'sid', 'siy', 'sin', 'sii', 'mse', 'aho', 'pai', 'pis', 'pcl', 'pas', 'tne', 'cma', 'tuz', 'vdi', 'tpe', 'asc', 'ffo', 'nmr', 'csy', 'nga', 'neq', 'kcr', 'hah']
#org_ids contains all oganisms supported by kaas 
values = {}

def init_values():

	return {'continue':mycontinue,'prog':prog,'uptype':uptype,'file':file,'peptide2':peptide2,\
'qname':qname,'mail':mail,'dbmode':'manual','org_list':org_list,'way':way,'mode':mode}

def write_res(response,outfile):

	outfile.write(response.read())
	
class Http_visiter(object):

	def __init__(self,url):
		
		self.url = url
		self.html = ''#created when using get or post function
		
	def get(self,url,values):
	
		data = urllib.urlencode(values)
		req = urr.Request(url+'?'+data)
		response = urr.urlopen(req)
		self.html = response.read()
		return self.html
		
	def post(self,url,values):
	
		register_openers()
		data, header = multipart_encode(values)
		req = urr.Request(url,data,header)
		response = urr.urlopen(req)
		self.html = response.read()
		return self.html
		
	def output(self,outfile):
	
		outfile.write(self.html)
	
	def refresh(self):
	
		self.__init__(self.url)
		
class KsRequest(Http_visiter):

	#kaas annot request
	
	def __init__(self,jobid='',status='',sttm='',edtm='',query_name="Query"):
		
		#sttm : start time
		#edtm : end time
		self.jobid = jobid
		self.status = status.strip()
		self.start_time = sttm
		self.end_time = edtm
		self.query_name = query_name
		self.res_header = 'http://www.genome.jp'
		self.res_url = ''
		
	def print_all(self):

		sys.stderr.write("\t".join([self.jobid,self.status,self.start_time,self.end_time,self.query_name])+"\n")
	
	def gen_url(self,par_url):
		
		#generate result url
		self.res_url = self.res_header+par_url
		
	def write_res(self,outfile):
	
		if self.res_url:
			values = ''
			outfile.write(self.get(self.res_url,values))
	
class Query_list(Http_visiter):
	
	#visit kegg query-list page convert to Query_list object
	
	def __init__(self,mail,url):

		super(Query_list,self).__init__(url)
		self.mail = mail
		self.values = [('mode','user'),('mail',self.mail)]
		self.get(self.url,self.values)
		self.data = self.parse(self.html)
		out = open('query_list.html','w')
		self.output(out)
		out.close()	
		
	def parse(self,html):
	
		requests = []
		q_soup = BeautifulSoup(html,"lxml")
		jobs = q_soup.find_all('tr')
		for child in jobs:#child eq a query
			tds = child.find_all('td')
			if tds[0].string:# == "ID"
				continue
			ks = KsRequest(tds[0].find('a').string,tds[2].string,tds[4].string,\
			tds[5].string,tds[1].string)
			if re.search('complete',ks.status,re.I):
				ks.gen_url(tds[3].find_all('a')[-1]["href"])
			requests.append(ks)
		return requests 
		
	def remove(self,id,mail):
	
		values = [('remove_id',id),('mail',mail),('mode','remove')]
		self.post(url,values)
		self.data = self.parse(self.html)
		return 
		
	def write_res(self,outfile,i=-1):
	
		self.data[i].write_res(outfile)
	
	def refresh(self):

		self.__init__(self.mail,self.url)
		
class Kaas_annot(Http_visiter):

	#parse compute response (socket._fileobject)
	url = 'http://www.genome.jp/kaas-bin/kaas_main'
	
	def __init__(self,values):
	
		#self.data = (status,jobid/description)
		#(status == 1 : job accepted,jobid)
		#(status == 0 : job refuesd,reason)
		super(Kaas_annot,self).__init__(Kaas_annot.url)
		self.values = values
		self.post(Kaas_annot.url,self.values)
		self.data = self.parse(self.html)
		
	def parse(self,html):
	
		c_soup = BeautifulSoup(html,"lxml")
		main = c_soup.find('div',id='main')
		mysub = main.find('p',class_='sub')
		if not mysub:
			status = 0
			jobid = main.find('p',class_='res').contents[0]
		elif re.search('Accepted',main.find('p',class_='sub').string,re.I):
			status = 1
			jobid = main.find('p',class_='res').string.split(":")[1].strip()
		elif re.search('Search program',main.find('p',class_='sub').string,re.I):
			sys.stderr.write("Error: one or more required information is missing\n")
			sys.exit()
		return (status,jobid)
		
	def jbinfo(self):
	
		#print job information
		if self.data[0]:
			sys.stderr.write("job accepted. job id : %s\n" % self.data[1])
		else:
			sys.stderr.write("job refuesd. reason : %s\n" % self.data[1])
	
	def get_status(self):
	
		return self.data[0]
		
	def get_jobid(self):

		if self.data[0]:
			return self.data[1]
		else:
			sys.stderr.write("job rejected.\n")
			
	def refresh(self):
	
		self.__init__(self.values)
		
if __name__ == '__main__':

	import sys, argparse, time
	
	parser = argparse.ArgumentParser(description='kaas-web api.It will be convenient to alter\
	the parameters [-p -m -g -w -n] in script')
	parser.add_argument('fastaFile',nargs='?',type=argparse.FileType('r'),\
			    default=sys.stdin,help='fastaFile [default: stdin]')
	parser.add_argument('-o','--outfile',nargs='?',type=argparse.FileType('w'),\
			    default=sys.stdout,help='outfile [default: stdout]')
	parser.add_argument('-p','--prog',nargs='?',help='b/x/z(BLAST/GHOSTX/GHOSTZ) [default:B]\
	 only when -t==p can you choose x/z',choices = ['b','x','z'])
	parser.add_argument('-m','--mail',nargs='?',help='e-mail address to receive result. \
	empty means you have revised in this script')
	parser.add_argument('-g','--org_list',nargs='?',help='hsa,dre,ath,etc.')
	parser.add_argument('-w','--way',nargs='?',help='b for BBH,s for SBH [default:b]',\
	choices = ['s','b'])
	parser.add_argument('-t','--seq_type',nargs='?',help='n for Nucleotide, p for protein [default:n]',\
	choices = ['n','p'])
	parser.add_argument('-r','--remove_jbid',nargs='?',help='type in job id')
	args = parser.parse_args(sys.argv[1:])
	
	if args.mail:
		mail = args.mail
	if mail == 'username@example.com':
		sys.stderr.write("please use your own e-mail.\n")
		sys.exit()
	if not re.match(r'\w+@\w+\.\w+',mail):
		sys.stderr.write("Error:bad mail format!\n")
		sys.exit()
		
	if args.remove_jbid:
		if not re.match(r'\d+$',args.remove_jbid):
			sys.stderr.write("Error:job id must be digit!\n")
			sys.exit()
		q = Query_list(mail,url)
		q.remove(args.remove_jbid,q.mail)
		sys.stderr.write("remove %s\n" % args.remove_jbid)
		sys.exit()
		
	if (args.prog == 'x' or args.prog == 'z') and args.seq_type != 'p':
		sys.stderr.write("Error:GHOSTX/GHOSTZ must be used with protein sequence\n\
		set -p x/z must with -t p\n")
		sys.exit()
	if args.prog == 'x':
		prog = 'GHOSTX'
	if args.prog == 'z':
		prog = 'GHOSTZ'
	if args.org_list:
		tmp = args.org_list.split(',')
		for each in tmp:
			if each.strip() not in org_ids:
				sys.stderr.write('Error: %s is not a valid organism\n'%each.strip())
				sys.exit()
		org_list = args.org_list
	if args.way == 's':
		way = 's'
	if args.seq_type == 'p':
		peptide2 = ''
		
	file = args.fastaFile
	values = init_values()
	
	def computing(q,outfile,cktm=0):
	
		#cktm checking time
		#cktm determine how long should it wait.
		assert isinstance(cktm,int)
		if cktm < 3:
			time.sleep(20)
		elif cktm < 6:
			time.sleep(180)
		else:
			time.sleep(600)
		cktm += 1
		q.refresh()
		mystatus = q.data[-1].status
		if re.search('complete',mystatus,re.I):
			sys.stderr.write("completed!!!\n")
			q.write_res(outfile)
		if re.search('failed',mystatus,re.I):
			sys.stderr.write("jobs failed.\n")
		if re.search('computing',mystatus,re.I):
			sys.stderr.write("computing...\n")
			computing(q,outfile,cktm)
			
	a = Kaas_annot(values)
	if not a.get_status(): #mission rejected
		q = Query_list(values['mail'],url)
		if re.match('computing',q.data[-1].status,re.I):
			sys.stderr.write("last work computing,please wait...\n")
		elif re.match('failed',q.data[-1].status,re.I):
			sys.stderr.write("last work failed.Remove it now. (job id:%s)\n" % q.data[-1].jobid)
			q.remove(q.data[-1].jobid,q.mail)
		else:
			sys.stderr.write("I don't know why this job is failed. \
			Please check your e-mail for more information\n")
	else: #mission accepted
		sys.stderr.write("job accepted.(job id : %s)\n" % a.get_jobid())
		sys.stderr.write("computing start\n")
		q = Query_list(values['mail'],url)
		computing(q,args.outfile)
