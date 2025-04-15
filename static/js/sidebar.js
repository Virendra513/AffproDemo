<script>
    document.addEventListener("DOMContentLoaded", function () {
        // Select all sidebar items
        const sidebarItems = document.querySelectorAll(".unique-item");
        const contentItems = document.querySelectorAll(".unique-content-item");

        sidebarItems.forEach(item => {
            item.addEventListener("click", function () {
                // Remove 'active' class from all items
                sidebarItems.forEach(i => i.classList.remove("active"));
                contentItems.forEach(c => c.classList.remove("active"));

                // Add 'active' class to the clicked item
                this.classList.add("active");

                // Get target content ID
                const targetId = this.getAttribute("data-target");
                document.getElementById(targetId).classList.add("active");
            });
        });
    });
</script>

