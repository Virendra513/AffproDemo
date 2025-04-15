document.addEventListener("DOMContentLoaded", function () {
    // Function to show content and hide others
    function showContent(sectionId) {
        document.querySelectorAll(".content").forEach(section => {
            section.style.display = "none";
        });
        document.getElementById(sectionId).style.display = "block";
    }

    // Event listener for sidebar links
    document.querySelectorAll(".sidebar-links a").forEach(link => {
        link.addEventListener("click", function (event) {
            event.preventDefault();
            let sectionId = this.getAttribute("onclick").match(/'([^']+)'/)[1];
            showContent(sectionId);
        });
    });

    // Sidebar toggle functionality
    const shrinkBtn = document.querySelector(".shrink-btn");
    const sidebar = document.querySelector("nav");
    const logoText = document.querySelector(".sidebar-top h3");
    
    shrinkBtn.addEventListener("click", function () {
        sidebar.classList.toggle("shrink");
        if (sidebar.classList.contains("shrink")) {
            logoText.classList.add("hide");
        } else {
            logoText.classList.remove("hide");
        }
    });

    // Restore last opened section (if needed)
    let lastSection = localStorage.getItem("lastSection") || "QClair-PLI";
    showContent(lastSection);

    document.querySelectorAll(".sidebar-links a").forEach(link => {
        link.addEventListener("click", function () {
            let sectionId = this.getAttribute("onclick").match(/'([^']+)'/)[1];
            localStorage.setItem("lastSection", sectionId);
        });
    });

    // Local storage handling for Predict section form
    const inputLigand = document.querySelector("#input1_pl");
    const inputProtein = document.querySelector("#input2_pl");
    const clearBtn = document.querySelector("#clearBtn");

    if (inputLigand && inputProtein && clearBtn) {
        inputLigand.value = localStorage.getItem("input1_pl") || "";
        inputProtein.value = localStorage.getItem("input2_pl") || "";
    
        document.querySelector("#projectForm").addEventListener("submit", function () {
            localStorage.setItem("input1_pl", inputLigand.value);
            localStorage.setItem("input2_pl", inputProtein.value);
        });
    
        clearBtn.addEventListener("click", function () {
            localStorage.removeItem("input1_pl");
            localStorage.removeItem("input2_pl");
            inputLigand.value = "";
            inputProtein.value = "";
        });
    }
});

