window.addEventListener("message", function(event) {
		try {
			var externalCall = JSON.parse(event.data);
			marvin.onReady(function() {
				marvin.sketcherInstance[externalCall.method].apply(marvin.sketcherInstance, externalCall.args);
			});
		} catch (e) {
			console.log(e);
		}
	}, false);


		// called when Marvin JS loaded
		function sketchOnLoad() {
			if(marvin.Sketch.isSupported()) {
				marvin.sketcherInstance = new marvin.Sketch("sketch");
				marvin.sketcherInstance.setServices(getDefaultServices());
			} else {
				alert("Cannot initiate sketcher. Current browser may not support HTML canvas or may run in Compatibility Mode.");
			}
		}
