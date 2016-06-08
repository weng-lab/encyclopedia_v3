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
        $("#b1").hide();
    });
    $("#b2").click(function(){
        $("#a2").css("visibility", "visible");
        $("#b2").hide();
    });
    $("#b3").click(function(){
        $("#a3").css("visibility", "visible");
        $("#b3").hide();
    });

    $("#gb1").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr18:55735848-56500847",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac":["hepg2_child_15_year"]
                 };
        processFormSubmitRet(event, "../enhancers/ucsc", fd);
    });

    $("#gb2").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr4:55090519-55159667",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac":["astrocyte_select"]
                 };
        processFormSubmitRet(event, "../enhancers/ucsc", fd);
    });

    $("#gb3").click(function(event){
        var fd = {"assembly":"hg19",
                  "loci":"chr16:52598938-52599438",
                  "assays":"BothDNaseAndH3K27ac",
                  "hg19BothDNaseAndH3K27ac":["mcf-7_adult_69_year"]
                 };
        processFormSubmitRet(event, "../enhancers/ucsc", fd);
    });

})
