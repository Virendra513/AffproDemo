document.addEventListener("DOMContentLoaded", function () {
    // Function to show content and hide others
    function showContent(sectionId) {
        document.querySelectorAll(".content").forEach(section => {
            section.style.display = "none";
        });
        document.getElementById(sectionId).style.display = "block";
    }

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

    // Sidebar links handling
    const sidebarLinks = document.querySelectorAll(".sidebar-links a");

    sidebarLinks.forEach(link => {
        link.addEventListener("click", function (event) {
            event.preventDefault(); // Prevent default link behavior
            let sectionId = this.getAttribute("onclick").match(/'([^']+)'/)[1];
            showContent(sectionId); // Show the corresponding section

            // Save the active tab in localStorage so we can restore it after form submission
            localStorage.setItem("lastSection", sectionId);
        });
    });

    // Highlight active tab and move active tab indicator
    let activeIndex;

    sidebarLinks.forEach((link, index) => {
        link.addEventListener("click", function () {
            // Remove the 'active' class from all sidebar links
            sidebarLinks.forEach(link => link.classList.remove("active"));
            // Add 'active' class to the clicked link
            this.classList.add("active");

            // Update activeIndex and move the active tab indicator
            activeIndex = index;
            moveActiveTab();
        });
    });

    // Function to move the active tab indicator
    function moveActiveTab() {
        const activeTab = document.querySelector(".active-tab");
        const shortcuts = document.querySelector(".sidebar-links h4");
        let topPosition = activeIndex * 58 + 2.5;

        if (activeIndex > 3) {
            topPosition += shortcuts.clientHeight;
        }

        activeTab.style.top = `${topPosition}px`;
    }

    // Local storage handling for Predict section form
    const inputLigand = document.querySelector("#input1_pl");
    const inputProtein = document.querySelector("#input2_pl");
    const clearBtn = document.querySelector("#clearBtn");

    if (inputLigand && inputProtein && clearBtn) {
        inputLigand.value = localStorage.getItem("input1_pl") || "";
        inputProtein.value = localStorage.getItem("input2_pl") || "";

        document.querySelector("#projectForm").addEventListener("submit", function (event) {
            localStorage.setItem("input1_pl", inputLigand.value);
            localStorage.setItem("input2_pl", inputProtein.value);

            // After submission, stay on the same tab (get last active section from localStorage)
            let lastSection = localStorage.getItem("lastSection") || "QClair-PLI";
            showContent(lastSection);

            // Prevent form from submitting (for demo purposes)
            event.preventDefault();
        });

        clearBtn.addEventListener("click", function () {
            localStorage.removeItem("input1_pl");
            localStorage.removeItem("input2_pl");
            inputLigand.value = "";
            inputProtein.value = "";
        });
    }

    // Restore last opened section (if needed)
    let lastSection = localStorage.getItem("lastSection") || "QClair-PLI";
    showContent(lastSection);
    sidebarLinks.forEach(link => {
        link.addEventListener("click", function () {
            let sectionId = this.getAttribute("onclick").match(/'([^']+)'/)[1];
            localStorage.setItem("lastSection", sectionId);
        });
    });
});
