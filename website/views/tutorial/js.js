function processFormSubmitRet(event, local_url, override_data){
    event.preventDefault();

    var formData = JSON.stringify(override_data)

    //console.log(formData);

    $.ajax({
        type: "POST",
        url: local_url,
        data: formData,
        dataType: "json",
        contentType : "application/json",
        async: false, // http://stackoverflow.com/a/20235765
        success: function(got){
            if("err" in got){
                $("#errMsg").text(got["err"]);
                $("#errBox").show()
                return true;
            }

            if("url" in got){
                return window.open(got["url"], '_blank');
            }
            var w = window.open();
            $(w.document.body).html(got["html"]);
        }
    });
}

$(document).ready(function(){
    $("#b1").click(function(){
        $("#a1").css("visibility", "visible");
    });

    $("#gb1").click(function(event){
        console.log(event);
        var fd = {"assembly":"hg19",
                  "loci":"chr18:55735848-56500847",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac":["hepg2_child_15_year"]
                 };
        processFormSubmitRet(event, "../enhancers/ucsc", fd);
    });

    $("#secondGalleryPane").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr17:40767000-40775000",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac": ["gm12878_select"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

    $("#thirdGalleryPane").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr6:42370000-42380000",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac": ["hela-s3_adult_31_year"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

    $("#fourthGalleryPane").click(function(event){
        var fd = {"assembly":"mm10",
                  "loci":"chr6:72214000-72220000",
                  "assays":"BothDNaseAndH3K27ac",
                  "mm10BothDNaseAndH3K27ac":["limb_embryonic_11_5_day"]
                 };
        processFormSubmitRet(event, "ucsc", fd);
    });

})
