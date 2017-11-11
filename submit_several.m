%% Submit a few jobs to the queue
for i=1:200
	mqsub('splineOFLIframePoints',{1001,-0.5,3,'pin',5,'xy'},'Partition','jila','Cores',4,'Name','ofli5xy','Time','4:00:00','Data','./Cluster','Memory',12)
end
