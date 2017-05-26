
; sets default environmental variables for NSTX-U

PRO set_genie_env

	MACHINE=getenv('MACHINE')
	IF MACHINE EQ '' THEN setenv, 'MACHINE=nstxu'	

	GENIE_PATH=getenv('GENIE_PATH')
	IF GENIE_PATH EQ '' THEN setenv, 'GENIE_PATH=/u/mreinke/GENIE/'

	ATOMIC_PHYSICS_PATH=getenv('ATOMIC_PHYSICS_PATH')
		IF ATOMIC_PHYSICS_PATH EQ '' THEN setenv, 'ATOMIC_PHYSICS_PATH=/p/nstxusr/nstx-users/mreinke/atomic_physics/'

	IMPSPEC_MDS_PATH=getenv('IMPSPEC_MDS_PATH')
	IF IMPSPEC_MDS_PATH EQ '' THEN $
	   setenv, 'IMPSPEC_MDS_PATH=\PASSIVESPEC::TOP.IMPSPEC'
	
	IMPSPEC_MDS_TREE=getenv('IMPSPEC_MDS_TREE')
	IF IMPSPEC_MDS_TREE EQ '' THEN $
	   setenv, 'IMPSPEC_MDS_TREE=passivespec'

	VIS_LINE_LIST=getenv('VIS_LINE_LIST')
	IF VIS_LINE_LIST EQ '' THEN $
	   setenv, 'VIS_LINE_LIST='+GENIE_PATH+'IMPSPEC/vis_line_list.tsv'
	
	VUV_LINE_LIST=getenv('VUV_LINE_LIST')
	IF VUV_LINE_LIST EQ '' THEN $
	   setenv, 'VUV_LINE_LIST='+GENIE_PATH+'IMPSPEC/vuv_line_list.tsv'
END

