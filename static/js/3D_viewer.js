
<script>
    var loadScriptAsync = function(uri) {
        return new Promise((resolve, reject) => {
            var tag = document.createElement('script');
            tag.src = uri;
            tag.async = true;
            tag.onload = resolve;
            tag.onerror = reject;
            document.head.appendChild(tag);
        });
    };

    var viewer_17378155808197615 = null;
    var warn = document.getElementById("3dmolwarning_17378155808197615");
    if (warn) {
        warn.parentNode.removeChild(warn);  // Remove warning message
    }

    // Load 3Dmol.js asynchronously
    loadScriptAsync('https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js')
        .then(function() {
            // Once 3Dmol.js is loaded, initialize the viewer
            viewer_17378155808197615 = $3Dmol.createViewer(document.getElementById("3dmolviewer_17378155808197615"), { backgroundColor: "white" });
            viewer_17378155808197615.zoomTo();
            viewer_17378155808197615.addModel("HETATM    1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00           C  \nHETATM    2  H1  UNL     1      -0.586  -0.190   0.923  1.00  0.00           H  \nHETATM    3  H2  UNL     1       0.460   1.008   0.056  1.00  0.00           H  \nHETATM    4  H3  UNL     1      -0.672  -0.052  -0.881  1.00  0.00           H  \nHETATM    5  H4  UNL     1       0.797  -0.766  -0.097  1.00  0.00           H  \nCONECT    1    2    3    4    5\nEND", "pdb");
            viewer_17378155808197615.setStyle({ "cartoon": { "color": "spectrum" } });
            viewer_17378155808197615.zoomTo();
            viewer_17378155808197615.render();
        })
        .catch(function(error) {
            console.error("Failed to load 3Dmol.js", error);
            alert("Error loading 3Dmol.js. Check the console for details.");
        });
</script>