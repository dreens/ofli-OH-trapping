%% Submit a few jobs to the queue
for i=1:500
	mqsub('splineOFLIframePoints',{1001,-0.5,3,'pin',5,'yz'},'Partition','jila','Cores',4,'Name','ofli5yz5','Time','4:00:00','Data','./Cluster','Memory',12)
end
for i=1:200
	mqsub('splineOFLIframePoints',{1001,-0.5,3,'pin',5,'yz'},'Partition','slow','Cores',4,'Name','ofli5yz5','Time','7-00:00:00','Data','./Cluster','Memory',8)
end
